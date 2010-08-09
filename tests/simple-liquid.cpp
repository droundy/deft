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
#include "Downhill.h"
#include "SteepestDescent.h"

const double kT = 1e-3; // room temperature in Hartree
const double ngas = 1e-5; // vapor density of water
const double mu = -kT*log(ngas);

int test_minimizer(const char *name, Minimizer *min, Grid *pot, int numiters, double fraccuracy=1e-3) {
  clock_t start = clock();
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  const double true_energy = -3.125598357241298;
  //const double gas_energy = -1.250000000000085e-11;

  *pot = +1e-4*((-10*pot->r2()).cwise().exp()) + 1.14*mu*VectorXd::Ones(pot->description().NxNyNz);

  for (int i=0;i<numiters && min->improve_energy(false);i++) {
    fflush(stdout);
  }

  min->print_info(-1);
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  printf("fractional energy error = %g\n", (min->energy() - true_energy)/fabs(true_energy));
  if (fabs((min->energy() - true_energy)/true_energy) > fraccuracy) {
    printf("Error in the energy is too big!\n");
    return 1;
  }
  return 0;
}

int main(int, char **argv) {
  Lattice lat(Cartesian(5,0,0), Cartesian(0,5,0), Cartesian(0,0,5));
  double resolution = 0.2;
  GridDescription gd(lat, resolution);
  Grid external_potential(gd);
  external_potential = 1000/ngas*(-0.2*external_potential.r2()).cwise().exp(); // repulsive bump
  Grid potential(gd);
  Functional attraction = GaussianPolynomial(gd, -1/ngas/ngas, 0.5, 2);
  Functional repulsion = GaussianPolynomial(gd, 1/ngas/ngas/ngas/ngas, 0.5, 4);
  Functional f0 = IdealGas(gd,kT) + ChemicalPotential(gd, mu) + ExternalPotential(external_potential)
    + attraction + repulsion;
  Functional f = compose(f0, EffectivePotentialToDensity(kT));

  Grid test_density(gd);
  test_density = EffectivePotentialToDensity(kT)(-1e-4*(-2*external_potential.r2()).cwise().exp() + mu*VectorXd::Ones(gd.NxNyNz));

  int retval = 0;

  potential = +1e-4*((-10*potential.r2()).cwise().exp()) + 1.14*mu*VectorXd::Ones(gd.NxNyNz);
  retval += f.run_finite_difference_test("simple liquid", potential);

  retval += attraction.run_finite_difference_test("quadratic", test_density);
  retval += repulsion.run_finite_difference_test("repulsive", test_density);

  Downhill downhill(f, &potential, 1e-3);
  potential.setZero();
  retval += test_minimizer("Downhill", &downhill, &potential, 5000, 0.995);

  PreconditionedDownhill pd(f, &potential, 1e-11);
  potential.setZero();
  retval += test_minimizer("PreconditionedDownhill", &pd, &potential, 300, 1e-14);

  SteepestDescent steepest(f, &potential, QuadraticLineMinimizer, 1e-3);
  potential.setZero();
  retval += test_minimizer("SteepestDescent", &steepest, &potential, 2000, 0.9);

  PreconditionedSteepestDescent psd(f, &potential, QuadraticLineMinimizer, 1e-11);
  potential.setZero();
  retval += test_minimizer("PreconditionedSteepestDescent", &psd, &potential, 200, 1e-10);

  
  potential = +1e-4*((-10*potential.r2()).cwise().exp()) + 1.14*mu*VectorXd::Ones(gd.NxNyNz);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
