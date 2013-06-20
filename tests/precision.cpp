// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Deft Authors
//
// Deft is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with deft; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include <stdio.h>
#include <time.h>
#include "OptimizedFunctionals.h"
#include "LineMinimizer.h"
#include "equation-of-state.h"

// Here we set up the lattice.
Lattice lat(Cartesian(5,0,0), Cartesian(0,5,0), Cartesian(0,0,5));
double resolution = 0.2;
GridDescription gd(lat, resolution);

Functional n = EffectivePotentialToDensity();
LiquidProperties prop = hughes_water_prop;

int test_minimizer(const char *name, Minimizer min, Functional f, Grid *pot,
                   double precision=1e-3, int maxiter = 100) {
  clock_t start = clock();
  printf("\n**********************************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s with precision %6g *\n", name, precision);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("**********************************\n\n");

  const double true_energy = -0.0003679537895448992;

  // First, let's reset the minimizer!
  min.minimize(f, pot->description(), pot);
  Minimizer foo = MaxIter(maxiter, Precision(precision, min));
  while (foo.improve_energy(false)) {
    //printf("Actual error is %g\n", foo.energy() - true_energy);
    fflush(stdout);
    //VectorXd pg = foo.pgrad();
    //ff.run_finite_difference_test("ff", *pot, &pg);
  }

  foo.print_info();
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  printf("energy error = %g\n", foo.energy() - true_energy);
  if (!(fabs(foo.energy() - true_energy) < precision)) { // double negatives handles NaNs correctly.
    printf("FAIL: Error in the energy is too big, should be %g! (%.16g vs %.16g)\n",
           precision, foo.energy(), true_energy);
    return 1;
  }
  return 0;
}

int main(int, char **argv) {
  int retval = 0;

  Grid potential(gd);
  Functional f = OfEffectivePotential(SaftFluid2(prop.lengthscale,
                                                 prop.epsilonAB, prop.kappaAB,
                                                 prop.epsilon_dispersion,
                                                 prop.lambda_dispersion, hughes_water_prop.length_scaling, 0));
  double mu = find_chemical_potential(f, prop.kT, prop.liquid_density);
  f = SaftFluid2(prop.lengthscale,
                 prop.epsilonAB, prop.kappaAB,
                 prop.epsilon_dispersion,
                 prop.lambda_dispersion, hughes_water_prop.length_scaling, mu);

  Grid external_potential(gd, 10*prop.kT*(-0.2*r2(gd)).cwise().exp()
                          - 5*prop.kT*(-0.05*r2(gd)).cwise().exp());
  Functional ff = OfEffectivePotential(f + ExternalPotential(external_potential));
  
  Minimizer downhill = Downhill(ff, gd, prop.kT, &potential, 1e-9);
  
  const double potential_value = -prop.kT*log(0.9*prop.liquid_density);
  potential.setZero();
  potential += potential_value*VectorXd::Ones(gd.NxNyNz);
  //retval += test_minimizer("Downhill", downhill, ff, &potential, 1e-13, 300);
  //retval += test_minimizer("Downhill", downhill, ff, &potential, 1e-10, 300);

  Minimizer pd = PreconditionedDownhill(ff, gd, prop.kT, &potential, 1e-11);
  potential.setZero();
  potential += potential_value*VectorXd::Ones(gd.NxNyNz);
  //retval += test_minimizer("PreconditionedDownhill", pd, ff, &potential, 1e-13, 300);

  Minimizer steepest = SteepestDescent(ff, gd, prop.kT, &potential, QuadraticLineMinimizer, 1e-3);
  potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  //retval += test_minimizer("SteepestDescent", steepest, ff, &potential, 1e-13, 20);
  //retval += test_minimizer("SteepestDescent", steepest, ff, &potential, 1e-3, 20);

  Minimizer psd = PreconditionedSteepestDescent(ff, gd, prop.kT, &potential, QuadraticLineMinimizer, 1e-11);
  potential.setZero();
  potential += potential_value*VectorXd::Ones(gd.NxNyNz);
  //retval += test_minimizer("PreconditionedSteepestDescent", psd, ff, &potential, 1e-13, 20);
  //retval += test_minimizer("PreconditionedSteepestDescent", psd, ff, &potential, 1e-5, 10);
  //retval += test_minimizer("PreconditionedSteepestDescent", psd, ff, &potential, 1e-2, 10);

  Minimizer cg = ConjugateGradient(ff, gd, prop.kT, &potential, QuadraticLineMinimizer, 1e-3);
  //potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  //retval += test_minimizer("ConjugateGradient", cg, ff, &potential, 1e-9, 200);
  potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  retval += test_minimizer("ConjugateGradient", cg, ff, &potential, 1e-5, 20);
  //potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  //retval += test_minimizer("ConjugateGradient", cg, ff, &potential, 1e-12, 400);
  potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  retval += test_minimizer("ConjugateGradient", cg, ff, &potential, 1e-3, 15);

  Minimizer pcg = PreconditionedConjugateGradient(ff, gd, prop.kT, &potential, QuadraticLineMinimizer, 1e-11);
  potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  retval += test_minimizer("PreconditionedConjugateGradient", pcg, ff, &potential, 1e-11, 20);
  //potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  //retval += test_minimizer("PreconditionedConjugateGradient", pcg, ff, &potential, 1e-9, 16);
  potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  retval += test_minimizer("PreconditionedConjugateGradient", pcg, ff, &potential, 1e-5, 11);
  potential = potential_value*VectorXd::Ones(gd.NxNyNz);
  retval += test_minimizer("PreconditionedConjugateGradient", pcg, ff, &potential, 1e-2, 3);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
