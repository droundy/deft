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
#include "Functionals.h"
#include "LineMinimizer.h"

const double my_kT = 1e-3; // room temperature in Hartree
const double ngas = 1.14e-7; // vapor density of water
const double nliquid = 4.9388942e-3; // density of liquid water
const double mu = -1.1*my_kT*log(ngas);
const double Veff_liquid = -my_kT*log(nliquid);

// Here we set up the lattice.
Lattice lat(Cartesian(5,0,0), Cartesian(0,5,0), Cartesian(0,0,5));
double resolution = 0.2;
GridDescription gd(lat, resolution);

// And the functional...
const double interaction_energy_scale = 0.01;
Functional attraction = GaussianPolynomial(-interaction_energy_scale/nliquid/nliquid/2, 0.5, 2);
Functional repulsion = GaussianPolynomial(interaction_energy_scale/nliquid/nliquid/nliquid/nliquid/4, 0.125, 4);
Functional f0 = ChemicalPotential(mu) + attraction + repulsion;
Functional n = EffectivePotentialToDensity();
Functional f = IdealGasOfVeff + f0(n);

Grid potential(gd);
Grid external_potential(gd, 1e-3/nliquid*(-0.2*r2(gd)).cwise().exp()); // repulsive bump

Functional ff = IdealGasOfVeff + (f0 + ExternalPotential(external_potential))(n);


int test_minimizer(const char *name, Minimizer min, Functional f, Grid *pot,
                   double precision=1e-3, int maxiter = 100) {
  clock_t start = clock();
  printf("\n**********************************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s with precision %6g *\n", name, precision);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("**********************************\n\n");

  const double true_energy = -0.2664377410430248;
  //const double gas_energy = -1.250000000000085e-11;

  *pot = +1e-4*((-10*r2(gd)).cwise().exp()) + 1.14*Veff_liquid*VectorXd::Ones(pot->description().NxNyNz);

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
  if (fabs(foo.energy() - true_energy) > precision) {
    printf("FAIL: Error in the energy is too big, should be %g! (%.16g vs %.16g)\n",
           precision, foo.energy(), true_energy);
    return 1;
  }
  return 0;
}

int main(int, char **argv) {
  int retval = 0;

  Minimizer downhill = Downhill(ff, gd, my_kT, &potential, 1e-11);
  
  potential.setZero();
  //retval += test_minimizer("Downhill", downhill, ff, &potential, 1e-13, 300);
  //retval += test_minimizer("Downhill", downhill, ff, &potential, 1e-10, 300);

  Minimizer pd = PreconditionedDownhill(ff, gd, my_kT, &potential, 1e-11);
  potential.setZero();
  //retval += test_minimizer("PreconditionedDownhill", pd, ff, &potential, 1e-13, 300);

  Minimizer steepest = SteepestDescent(ff, gd, my_kT, &potential, QuadraticLineMinimizer, 1e-3);
  potential.setZero();
  //retval += test_minimizer("SteepestDescent", steepest, ff, &potential, 1e-13, 20);
  //retval += test_minimizer("SteepestDescent", steepest, ff, &potential, 1e-3, 20);

  Minimizer psd = PreconditionedSteepestDescent(ff, gd, my_kT, &potential, QuadraticLineMinimizer, 1e-11);
  potential.setZero();
  //retval += test_minimizer("PreconditionedSteepestDescent", psd, ff, &potential, 1e-13, 20);
  //retval += test_minimizer("PreconditionedSteepestDescent", psd, ff, &potential, 1e-5, 10);
  //retval += test_minimizer("PreconditionedSteepestDescent", psd, ff, &potential, 1e-2, 10);

  Minimizer cg = ConjugateGradient(ff, gd, my_kT, &potential, QuadraticLineMinimizer, 1e-3);
  potential.setZero();
  retval += test_minimizer("ConjugateGradient", cg, ff, &potential, 1e-12, 90);
  retval += test_minimizer("ConjugateGradient", cg, ff, &potential, 1e-3, 20);

  Minimizer pcg = PreconditionedConjugateGradient(ff, gd, my_kT, &potential, QuadraticLineMinimizer, 1e-11);
  potential.setZero();
  // Oddly enough with PreconditionedConjugateGradient we can't get
  // better than 1e-11 precision, and even that requires setting the
  // Precision requirement extra high.
  retval += test_minimizer("PreconditionedConjugateGradient", pcg, ff, &potential, 1e-11, 60);
  retval += test_minimizer("PreconditionedConjugateGradient", pcg, ff, &potential, 1e-5, 10);
  retval += test_minimizer("PreconditionedConjugateGradient", pcg, ff, &potential, 1e-2, 10);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
