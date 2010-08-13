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
const double ngas = 1.14e-7; // vapor density of water
const double nliquid = 4.9388942e-3; // density of liquid water
const double mu = -kT*log(ngas);
const double Veff_liquid = -kT*log(nliquid);

// Here we set up the lattice.
Lattice lat(Cartesian(5,0,0), Cartesian(0,5,0), Cartesian(0,0,5));
double resolution = 0.2;
GridDescription gd(lat, resolution);

// And the functional...
const double interaction_energy_scale = 0.01;
Functional attraction = GaussianPolynomial(-interaction_energy_scale/nliquid/nliquid/2, 0.5, 2);
Functional repulsion = GaussianPolynomial(interaction_energy_scale/nliquid/nliquid/nliquid/nliquid/4, 0.125, 4);
Functional f0 = IdealGas(kT) + ChemicalPotential(mu) + attraction + repulsion;
FieldFunctional n = EffectivePotentialToDensity(kT);
Functional f = f0(n);

Grid potential(gd);
Grid external_potential(gd, 1e-3/nliquid*(-0.2*potential.r2()).cwise().exp()); // repulsive bump

Functional ff = (f0 + ExternalPotential(external_potential))(n);


int test_minimizer(const char *name, Minimizer min, Grid *pot, double fraccuracy=1e-3) {
  clock_t start = clock();
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  const double true_energy = -0.2639209612386242;
  //const double gas_energy = -1.250000000000085e-11;

  *pot = +1e-4*((-10*pot->r2()).cwise().exp()) + 1.14*Veff_liquid*VectorXd::Ones(pot->description().NxNyNz);

  while (min.improve_energy(false)) fflush(stdout);

  min.print_info();
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  printf("fractional energy error = %g\n", (min.energy() - true_energy)/fabs(true_energy));
  if (fabs((min.energy() - true_energy)/true_energy) > fraccuracy) {
    printf("Error in the energy is too big!\n");
    return 1;
  }
  return 0;
}

int main(int, char **argv) {
  int retval = 0;

  {
    Grid test_density(gd, EffectivePotentialToDensity(kT)(gd, -1e-4*(-2*external_potential.r2()).cwise().exp()
                                                          + mu*VectorXd::Ones(gd.NxNyNz)));
    potential = +1e-4*((-10*potential.r2()).cwise().exp()) + 1.14*Veff_liquid*VectorXd::Ones(gd.NxNyNz);
    retval += f.run_finite_difference_test("simple liquid", potential);
    
    retval += attraction.run_finite_difference_test("quadratic", test_density);
    retval += repulsion.run_finite_difference_test("repulsive", test_density);
  }

  Minimizer downhill = MaxIter(300, Downhill(ff, gd, &potential, 1e-11));
  potential.setZero();
  retval += test_minimizer("Downhill", downhill, &potential, 1e-13);

  Minimizer pd = MaxIter(300, PreconditionedDownhill(ff, gd, &potential, 1e-11));
  potential.setZero();
  retval += test_minimizer("PreconditionedDownhill", pd, &potential, 1e-13);

  Minimizer steepest = MaxIter(200, SteepestDescent(ff, gd, &potential, QuadraticLineMinimizer, 1e-3));
  potential.setZero();
  retval += test_minimizer("SteepestDescent", steepest, &potential, 1e-13);

  Minimizer psd = MaxIter(200, PreconditionedSteepestDescent(ff, gd, &potential, QuadraticLineMinimizer, 1e-11));
  potential.setZero();
  retval += test_minimizer("PreconditionedSteepestDescent", psd, &potential, 1e-13);

  
  potential = +1e-4*((-10*potential.r2()).cwise().exp()) + 1.14*mu*VectorXd::Ones(gd.NxNyNz);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
