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
#include "equation-of-state.h"
#include "LineMinimizer.h"

#include "utilities.h"

const double R = 2.7;
const double eta_one = 3.0/(4*M_PI*R*R*R);
const double diameter_cubed = 1/(8*R*R*R);
const double nliquid = 0.4257*eta_one;
const double mu = find_chemical_potential(OfEffectivePotential(HardSpheres(R) + IdealGas()),
                                          water_prop.kT, nliquid);

// Here we set up the lattice.
const double zmax = 80;
const double width = 0.0001;
Lattice lat(Cartesian(width,0,0), Cartesian(0,width,0), Cartesian(0,0,zmax));
GridDescription gd(lat, 0.1);

// And the functional...
Functional f0wb = HardSpheresWB(R);
Functional f0rf = HardSpheresRF(R);
Functional f = OfEffectivePotential(HardSphereGas(R, mu));
//Functional fslow = OfEffectivePotential(HardSpheres(R) + IdealGas() + ChemicalPotential(mu));
Functional fslow = HardSpheres(R)(EffectivePotentialToDensity())
  + IdealGas()(EffectivePotentialToDensity()) + ChemicalPotential(mu)(EffectivePotentialToDensity());
//Functional ig = IdealGas()(EffectivePotentialToDensity());
Functional ig = OfEffectivePotential(IdealGas());
Functional igveff = IdealGasOfVeff();

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

  potential = external_potential + 0.005*VectorXd::Ones(gd.NxNyNz);

  for (int i=0;i<numiters && min.improve_energy(false);i++) {
    fflush(stdout);
  }
  min.print_info();
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);

  const double true_energy = -0.0028439564619177;
  const double true_N = 0.376241423570245;

  int retval = 0;
  double energy = min.energy()/width/width;
  printf("Energy is %.15g\n", energy);
  if (energy < true_energy) {
    printf("FAIL: Energy is less than the true energy by %g!\n", true_energy - energy);
    retval++;
  }
  if (fabs(energy/true_energy - 1) > fraccuracy) {
    printf("FAIL: Error in the energy is too big: %g\n", (energy - true_energy)/fabs(true_energy));
    retval++;
  }

  double N = 0;
  {
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    for (int i=0;i<gd.NxNyNz;i++) N += density[i]*gd.dvolume;
  }
  N = N/width/width;
  printf("N is %.15g\n", N);
  if (fabs(N/true_N - 1) > 10*fraccuracy) {
    printf("FAIL: Error in N is too big: %g\n", N/true_N - 1);
    retval++;
  }
  return retval;
}

double walls(Cartesian r) {
  const double z = r.dot(Cartesian(0,0,1));
  if (fabs(z) < 2*R) {
    return 0.1;
  } else {
    return 0;
  }
}

double notinwall(Cartesian r) {
  const double z = r.dot(Cartesian(0,0,1));
  if (fabs(z) < 2*R) {
    return 0;
  } else {
    return 1;
  }
}

int main(int, char **argv) {
  printf("mu is %g\n", mu);
  external_potential.Set(walls);
  Grid constraint(gd);
  constraint.Set(notinwall);
  //Functional f1 = f0 + ExternalPotential(external_potential);
  ff = constrain(constraint, f);

  int retval = 0;

  {
    potential = external_potential + 0.005*VectorXd::Ones(gd.NxNyNz);
    retval += OfEffectivePotential(f0wb).run_finite_difference_test("white bear functional",
                                                                    water_prop.kT, potential);
    retval += OfEffectivePotential(f0rf).run_finite_difference_test("rosenfeld functional",
                                                                    water_prop.kT, potential);
    retval += ff.run_finite_difference_test("full functional", water_prop.kT, potential);
    if (retval > 0) {
      printf("quitting early with a failure in %d finite difference tests.\n", retval);
      return retval;
    }
  }

  {
    Minimizer pd = Precision(0, ConjugateGradient(ff, gd, water_prop.kT, &potential,
                                                  QuadraticLineMinimizer));
    retval += test_minimizer("ConjugateGradient", pd, 120, 1e-7);
  }
  {
    Minimizer pd =
      Precision(0, PreconditionedConjugateGradient(ff, gd, water_prop.kT, &potential, QuadraticLineMinimizer));
    retval += test_minimizer("PreconditionedConjugateGradient", pd, 95, 1e-7);

    //potential = external_potential + mu*VectorXd::Ones(gd.NxNyNz);
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    //density.epsNativeSlice("PreconditionedConjugateGradient.eps", Cartesian(0,0,zmax), Cartesian(1,0,0),
    //                       Cartesian(0,0,0));
    density.epsNative1d("hard-wall-plot.eps", Cartesian(0,0,0), Cartesian(0,0,zmax), diameter_cubed, R, "Y axis: n*8*R^3, x axis: R");
    potential.epsNative1d("hard-wall-potential.eps", Cartesian(0,0,0), Cartesian(0,0,zmax), 1, R);
    Grid grad(gd), pgrad(gd);
    grad.setZero();
    pgrad.setZero();
    ff.integralgrad(water_prop.kT, potential, &grad, &pgrad);
    grad.epsNative1d("hard-wall-grad.eps", Cartesian(0,0,0), Cartesian(0,0,zmax), 1, R);
    pgrad.epsNative1d("hard-wall-pgrad.eps", Cartesian(0,0,0), Cartesian(0,0,zmax), 1, R);
    Grid(gd, StepConvolve(R)(water_prop.kT, density)).epsNative1d("n3.eps",
                                                                  Cartesian(0,0,0), Cartesian(0,0,zmax), 1, R);
 
    retval += constrain(constraint, f0wb).run_finite_difference_test("white bear functional",
                                                                     water_prop.kT, density, &grad);
    retval += constrain(constraint, f0rf).run_finite_difference_test("rosenfeld functional",
                                                                     water_prop.kT, density, &grad);
  }

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
