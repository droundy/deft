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

const double kT = 1e-3; // room temperature in Hartree
const double R = 3.0;
const double ngas = 0.005; // approximate density of interest
const double mu = -kT*log(ngas);

// Here we set up the lattice.
Lattice lat(Cartesian(10,0,0), Cartesian(0,10,0), Cartesian(0,0,10));
double resolution = 0.5;
GridDescription gd(lat, resolution);

// And the functional...
const double interaction_energy_scale = 0.01;
Functional f00 = HardSpheres(R, kT);
Functional f0 = IdealGas(kT) + f00 + ChemicalPotential(mu);
FieldFunctional n = EffectivePotentialToDensity(kT);
Functional f = f0(n);

Grid potential(gd);
Grid external_potential(gd, 1e-3/ngas*(-0.2*potential.r2()).cwise().exp()); // repulsive bump

Functional ff = (f0 + ExternalPotential(external_potential))(n);


int test_minimizer(const char *name, Minimizer min, Grid *pot, double fraccuracy=1e-3) {
  clock_t start = clock();
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  const double true_energy = -0.2639034579484159;
  //const double gas_energy = -1.250000000000085e-11;

  *pot = +1e-4*((-10*pot->r2()).cwise().exp()) + 1.14*mu*VectorXd::Ones(pot->description().NxNyNz);

  while (min.improve_energy(true)) fflush(stdout);

  min.print_info();
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  printf("fractional energy error = %g\n", (min.energy() - true_energy)/fabs(true_energy));
  if (fabs((min.energy() - true_energy)/true_energy) > fraccuracy) {
    printf("Error in the energy is too big!\n");
    return 1;
  }
  if (min.energy() < true_energy) {
    printf("Sign of error is wrong!!!\n");
    return 1;
  }
  return 0;
}

int main(int, char **argv) {
  int retval = 0;

  {
    Grid test_density(gd, EffectivePotentialToDensity(kT)(gd, -0.04*(-2*external_potential.r2()).cwise().exp()
                                                          + mu*VectorXd::Ones(gd.NxNyNz)));
    potential = +1e-4*((-10*potential.r2()).cwise().exp()) + 1.14*mu*VectorXd::Ones(gd.NxNyNz);
    //retval += f00.run_finite_difference_test("hard spheres with no ideal gas", test_density);
    //retval += f0.run_finite_difference_test("hard spheres straight", test_density);
    const double expected_energy = 0.002792960896652312;
    printf("hard sphere energy is %.16g\n", f(potential));
    if (fabs(f(potential)/expected_energy - 1) > 1e-14) {
      printf("Error in hard sphere energy (of fixed potential) is too big: %g (from %.16g)\n",
             f(potential)/expected_energy - 1, f(potential));
      retval++;
    }
    retval += f.run_finite_difference_test("hard spheres", potential);
  }

  /*
  Minimizer downhill = MaxIter(300, Downhill(ff, gd, &potential));
  potential.setZero();
  retval += test_minimizer("Downhill", downhill, &potential, 1e-9);

  Minimizer pd = MaxIter(300, PreconditionedDownhill(ff, gd, &potential));
  potential.setZero();
  retval += test_minimizer("PreconditionedDownhill", pd, &potential, 1e-9);

  Minimizer steepest = MaxIter(100, SteepestDescent(ff, gd, &potential, QuadraticLineMinimizer));
  potential.setZero();
  retval += test_minimizer("SteepestDescent", steepest, &potential, 1e-9);

  Minimizer psd = MaxIter(100, PreconditionedSteepestDescent(ff, gd, &potential, QuadraticLineMinimizer));
  potential.setZero();
  retval += test_minimizer("PreconditionedSteepestDescent", psd, &potential, 1e-9);


  Minimizer cg = MaxIter(5, ConjugateGradient(ff, gd, &potential, QuadraticLineMinimizer));
  potential.setZero();
  retval += test_minimizer("ConjugateGradient", cg, &potential, 1e-15); // I used this to verify...

  Grid density(gd, EffectivePotentialToDensity(kT)(gd, potential));
  density.epsNativeSlice("hardspheres.eps", Cartesian(10,0,0), Cartesian(0,10,0),
                         Cartesian(0,0,0));
  printf("Energy is really %g\n", (f0 + ExternalPotential(external_potential))(density));
  Grid n3(gd, StepConvolve(R)(density));
  n3.epsNativeSlice("n3.eps", Cartesian(10,0,0), Cartesian(0,10,0), Cartesian(0,0,0));

  //Minimizer pcg = MaxIter(100, PreconditionedConjugateGradient(ff, gd, &potential, QuadraticLineMinimizer));
  //potential.setZero();
  //retval += test_minimizer("PreconditionedConjugateGradient", pcg, &potential, 1e-11);
  */

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
