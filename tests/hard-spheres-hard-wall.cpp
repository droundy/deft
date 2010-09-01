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
#include "utilities.h"

const double kT = water_prop.kT; // room temperature in Hartree
const double R = 2.7;
const double nliquid = water_prop.liquid_density; // density of liquid water
const double mu = -kT*log(nliquid);

// Here we set up the lattice.
const double zmax = 40;
Lattice lat(Cartesian(0.01,0,0), Cartesian(0,0.01,0), Cartesian(0,0,zmax));
GridDescription gd(lat, 0.1);

// And the functional...
Functional f00 = HardSpheres(R, kT);
Functional f0 = IdealGas(kT) + f00 + ChemicalPotential(mu);
FieldFunctional n = EffectivePotentialToDensity(kT);
Functional f = f0(n);

Grid external_potential(gd);
Grid potential(gd);
Functional ff;

int test_minimizer(const char *name, Minimizer min, int numiters, double fraccuracy=1e-3) {
  clock_t start = clock();
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  potential = external_potential + mu*VectorXd::Ones(gd.NxNyNz);

  for (int i=0;i<numiters && min.improve_energy(true);i++) {
    fflush(stdout);
  }
  min.print_info();
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);

  const double true_energy = -9.26646560514309e-09;
  const double true_N = 5.97351547018021e-06;

  int retval = 0;
  printf("Energy is %.15g\n", min.energy());
  if (min.energy() < true_energy) {
    printf("FAIL: Energy is less than the true energy by %g!\n", true_energy - min.energy());
    retval++;
  }
  if (fabs(min.energy()/true_energy - 1) > fraccuracy) {
    printf("FAIL: Error in the energy is too big: %g\n", (min.energy() - true_energy)/fabs(true_energy));
    retval++;
  }

  double N = 0;
  {
    Grid density(gd, EffectivePotentialToDensity(kT)(gd, potential));
    for (int i=0;i<gd.NxNyNz;i++) N += density[i]*gd.dvolume;
  }
  printf("N is %.15g\n", N);
  if (fabs(N/true_N - 1) > fraccuracy) {
    printf("FAIL: Error in N is too big: %g\n", N/true_N - 1);
    retval++;
  }
  return retval;
}

double walls(Cartesian r) {
  const double z = r.dot(Cartesian(0,0,1));
  if (fabs(z) < 1.5*R) {
    return 1e9;
  } else {
    return 0;
  }
}

double notinwall(Cartesian r) {
  const double z = r.dot(Cartesian(0,0,1));
  if (fabs(z) < 1.5*R) {
    return 0;
  } else {
    return 1;
  }
}

int main(int, char **argv) {
  external_potential.Set(walls);
  Grid constraint(gd);
  constraint.Set(notinwall);
  //Functional f1 = f0 + ExternalPotential(external_potential);
  Functional f1 = constrain(constraint, f0);
  ff = f1(n);

  int retval = 0;
  {
    Minimizer pd = Precision(0, PreconditionedConjugateGradient(ff, gd, &potential, QuadraticLineMinimizer));
    retval += test_minimizer("PreconditionedConjugateGradient", pd, 20, 1e-5);

    Grid density(gd, EffectivePotentialToDensity(kT)(gd, potential));
    density.epsNativeSlice("PreconditionedConjugateGradient.eps", Cartesian(0,0,zmax), Cartesian(1,0,0),
                           Cartesian(0,0,0));
    density.epsNative1d("hard-wall-plot.eps", Cartesian(0,0,0), Cartesian(0,0,zmax));

    //retval += f0.run_finite_difference_test("f0", density);
    //retval += f00.run_finite_difference_test("f00", density);
  }

  //Minimizer psd = PreconditionedSteepestDescent(ff, gd, &potential, QuadraticLineMinimizer, 1e-4);
  //retval += test_minimizer("PreconditionedSteepestDescent", psd, 20, 1e-4);

  //Minimizer cg = ConjugateGradient(ff, gd, &potential, QuadraticLineMinimizer);
  //retval += test_minimizer("ConjugateGradient", cg, 3000, 3e-5);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
