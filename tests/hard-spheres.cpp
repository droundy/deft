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
const double R = 2.7;
const double ngas = 0.005; // approximate density of interest
const double mu = -kT*log(ngas);

// Here we set up the lattice.
Lattice lat(Cartesian(10,0,0), Cartesian(0,10,0), Cartesian(0,0,10));
double resolution = 0.5;
GridDescription gd(lat, resolution);

// And the functional...
Functional f00 = integrate(HardSpheres(R, kT));
Functional f0 = integrate(IdealGas(kT) + HardSpheres(R, kT) + ChemicalPotential(mu));
FieldFunctional n = EffectivePotentialToDensity(kT);
Functional f = f0(n);

Grid potential(gd);
Grid external_potential(gd, 1e-3/ngas*(-0.2*potential.r2()).cwise().exp()); // repulsive bump

Functional ff = (f0 + integrate(ExternalPotential(external_potential)))(n);


int test_minimizer(const char *name, Minimizer min, Grid *pot, double accuracy=1e-3) {
  clock_t start = clock();
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  const double true_energy = -0.0013480883207912;

  *pot = +1e-4*((-10*pot->r2()).cwise().exp()) + 1.14*mu*VectorXd::Ones(pot->description().NxNyNz);

  while (min.improve_energy(true)) {
    printf("Actual error is %g\n", min.energy() - true_energy);
    fflush(stdout);
  }

  min.print_info();
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  printf("energy error = %g\n", min.energy() - true_energy);
  if (min.energy() - true_energy > accuracy) {
    printf("FAIL: Error in the energy is too big!\n");
    return 1;
  }
  if (min.energy() < true_energy) {
    printf("FAIL: Sign of error is wrong!!!\n");
    return 1;
  }
  return 0;
}

int main(int, char **argv) {
  int retval = 0;

  {
    potential = +1e-4*((-10*potential.r2()).cwise().exp()) + 1.14*mu*VectorXd::Ones(gd.NxNyNz);
    Grid test_density(gd, EffectivePotentialToDensity(kT)(gd, potential));
    //Grid test_density(gd, EffectivePotentialToDensity(kT)(gd, -0.04*(-2*external_potential.r2()).cwise().exp()
    //                                                      + mu*VectorXd::Ones(gd.NxNyNz)));
    //retval += f00.run_finite_difference_test("hard spheres with no ideal gas", test_density);
    //retval += f0.run_finite_difference_test("hard spheres straight", test_density);
    const double expected_energy = -0.002709673141096069;
    printf("hard sphere energy is %.16g\n", f(potential));
    if (fabs(f(potential)/expected_energy - 1) > 1e-13) {
      printf("FAIL: Error in hard sphere energy (of fixed potential) is too big %g (from %.16g)\n",
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
  */

  Minimizer steepest = Precision(1e-5, MaxIter(20, SteepestDescent(ff, gd, &potential, QuadraticLineMinimizer)));
  potential.setZero();
  retval += test_minimizer("SteepestDescent", steepest, &potential, 3e-5);

  Minimizer psd = Precision(1e-6, MaxIter(15, PreconditionedSteepestDescent(ff, gd, &potential, QuadraticLineMinimizer)));
  potential.setZero();
  retval += test_minimizer("PreconditionedSteepestDescent", psd, &potential, 3e-6);

  Minimizer cg = Precision(1e-6, MaxIter(25, ConjugateGradient(ff, gd, &potential, QuadraticLineMinimizer)));
  potential.setZero();
  retval += test_minimizer("ConjugateGradient", cg, &potential, 3e-6);

  Grid density(gd, EffectivePotentialToDensity(kT)(gd, potential));
  density.epsNativeSlice("hardspheres.eps", Cartesian(10,0,0), Cartesian(0,10,0), Cartesian(0,0,0));

  Minimizer pcg = Precision(1e-6, MaxIter(15, PreconditionedConjugateGradient(ff, gd, &potential, QuadraticLineMinimizer)));
  potential.setZero();
  retval += test_minimizer("PreconditionedConjugateGradient", pcg, &potential, 3e-6);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
