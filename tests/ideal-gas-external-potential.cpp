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
#include "Functionals.h"
#include "LineMinimizer.h"

const double kT = 1e-3; // room temperature in Hartree
const double ngas = 1e-5; // vapor density of water
const double mu = -kT*log(ngas);

int test_minimizer(const char *name, Minimizer min, Grid *pot, Grid expected_density, int numiters) {
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  const double true_energy = -5.597856610022806e-11;

  *pot = +1e-4*((-10*pot->r2()).cwise().exp()) + 1.04*mu*VectorXd::Ones(pot->description().NxNyNz);
  for (int i=0;i<numiters && min.improve_energy(false);i++) {
    fflush(stdout);
  }
  min.print_info();
  Grid density(pot->description(), EffectivePotentialToDensity(kT)(*pot));
  double err2 = 0;
  for (int i=0;i<pot->description().NxNyNz;i++) {
    err2 += (density[i]-expected_density[i])*(density[i]-expected_density[i]);
  }
  err2 /= pot->description().NxNyNz;
  printf("rms error = %g\n", sqrt(err2));
  printf("fractional energy error = %g\n", (min.energy() - true_energy)/fabs(true_energy));
  if (fabs((min.energy() - true_energy)/true_energy) > 5e-3) {
    printf("Error in the energy is too big!\n");
    return 1;
  }
  for (int i=0;i<pot->description().NxNyNz;i++) {
    if (fabs(density[i] - expected_density[i]) > fabs(1e-3*ngas)) {
      printf("Oh no, the error in density is %g out of %g at %d (%g)!\n",
             density[i] - expected_density[i], expected_density[i], i,
             fabs(density[i] - expected_density[i])/ngas);
      return 1;
    }
  }
  return 0;
}

int main(int, char **argv) {
  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  int resolution = 5;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid external_potential(gd);
  external_potential = 1e-1*external_potential.r2(); // a harmonic trap...
  Grid potential(gd);
  //potential.epsNativeSlice("potential.eps", Cartesian(1,0,0),
  //                         Cartesian(0,1,0), Cartesian(0,0,0));
  Functional ig_and_mu = IdealGas(kT) + ChemicalPotential(mu) + ExternalPotential(external_potential);
  Functional f = compose(ig_and_mu, EffectivePotentialToDensity(kT));

  Grid expected_density(gd, EffectivePotentialToDensity(kT)(gd, external_potential + mu*VectorXd::Ones(gd.NxNyNz)));

  int retval = 0;

  Minimizer downhill = Downhill(f, gd, &potential);
  potential.setZero();
  retval += test_minimizer("Downhill", downhill, &potential, expected_density, 3000);

  Minimizer pd = PreconditionedDownhill(f, gd, &potential);
  potential.setZero();
  retval += test_minimizer("PreconditionedDownhill", pd, &potential, expected_density, 15);

  Minimizer steepest = SteepestDescent(f, gd, &potential, QuadraticLineMinimizer, 1.0);
  potential.setZero();
  retval += test_minimizer("SteepestDescent", steepest, &potential, expected_density, 2000);

  Minimizer psd = PreconditionedSteepestDescent(f, gd, &potential, QuadraticLineMinimizer, 1.0);
  potential.setZero();
  retval += test_minimizer("PreconditionedSteepestDescent", psd, &potential, expected_density, 10);

  if (retval == 0) {
    printf("%s passes!\n", argv[0]);
  } else {
    printf("%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
