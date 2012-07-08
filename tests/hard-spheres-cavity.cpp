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
#include "utilities.h"

const double R = 2.7;
const double eta_one = 3.0/(4*M_PI*R*R*R);
const double diameter_cubed = 1/(8*R*R*R);
const double nliquid = 0.324*eta_one;
Functional n = EffectivePotentialToDensity();
const double mu = find_chemical_potential(HardSpheres(R)(n) + IdealGasOfVeff(), water_prop.kT, nliquid);

// Here we set up the lattice.
const double rcav = R+R; // 11.8*R+R;
const double rmax = rcav*2;
Lattice lat(Cartesian(0,rmax,rmax), Cartesian(rmax,0,rmax), Cartesian(rmax,rmax,0));
//Lattice lat(Cartesian(1.4*rmax,0,0), Cartesian(0,1.4*rmax,0), Cartesian(0,0,1.4*rmax));
GridDescription gd(lat, 0.5);

Grid external_potential(gd);
Grid potential(gd);

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

  const double true_energy = -0.0375954069264892;

  int retval = 0;
  double energy = min.energy();
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
    Grid density(gd, n(water_prop.kT, gd, potential));
    for (int i=0;i<gd.NxNyNz;i++) N += density[i]*gd.dvolume;
  }
  printf("N is %.15g\n", N);
  //if (fabs(N/true_N - 1) > 10*fraccuracy) {
  //  printf("FAIL: Error in N is too big: %g\n", N/true_N - 1);
  //  retval++;
  //}
  return retval;
}

double notincavity(Cartesian r) {
  const double rad2 = r.dot(r);
  if (rad2 < rcav*rcav) {
    return 0;
  } else {
    return 1;
  }
}

double incavity(Cartesian r) {
  return 1 - notincavity(r);
}

void this_took() {
  static clock_t start = 0;
  clock_t end = clock();
  if (start)
    printf("    This took %g seconds.\n", (end - double(start))/CLOCKS_PER_SEC);
  start = end;
}

int main(int, char **argv) {
  Grid constraint(gd);
  this_took();
  printf("I am about to set the constraint...\n"); 
  constraint.Set(notincavity);
  this_took();
  printf("I have set the constraint...\n");
  // The functionals are...
  Functional n = constrain(constraint, EffectivePotentialToDensity());
  Functional f0wb = HardSpheresWB(R);
  Functional f0rf = HardSpheresRFFast(R);
  Functional ff = constrain(constraint, OfEffectivePotential(HardSphereGas(R, mu)));

  printf("I am about to set the initial cavity...\n"); 
  external_potential.Set(incavity);
  external_potential *= 0.1;
  printf("I have set the initial cavity...\n");
  this_took();
  external_potential.epsNativeSlice("external.eps", Cartesian(2*rmax,0,0), Cartesian(0,2*rmax,0), Cartesian(-rmax,-rmax,0));
  printf("I have output a native slice of the external potential...\n");
  this_took();
  external_potential.epsRadial1d("external-radial.eps", 0, rmax, 1, R, "Good fun!");
  printf("I have output a radial slice of the external potential...\n");
  this_took();

  int retval = 0;

  {
    GridDescription gd(lat, 0.5);
    Grid small_potential(gd);
    Functional n = EffectivePotentialToDensity();
    small_potential.Set(incavity);
    small_potential *= 1e9;
    small_potential = small_potential + 0.005*VectorXd::Ones(gd.NxNyNz);
    retval += f0wb(n).run_finite_difference_test("white bear functional", water_prop.kT, small_potential);
    retval += f0rf(n).run_finite_difference_test("rosenfeld functional", water_prop.kT, small_potential);
    printf("Done with both finite difference tests!\n");
  }

  {
    reset_peak_memory();
    {
      Minimizer pd = Precision(0, ConjugateGradient(ff, gd, water_prop.kT, &potential, QuadraticLineMinimizer));
      retval += test_minimizer("ConjugateGradient", pd, 100, 1e-5);

      double peak = peak_memory()/1024.0/1024;
      reset_peak_memory();
      printf("Peak memory use was %g\n", peak);
      if (peak > 7 || peak < 6) {
        printf("FAIL: Peak memory use was not as expected in ConjugateGradient...\n");
        retval++;
      }
    }
    {
      Minimizer pd = Precision(0, PreconditionedConjugateGradient(ff, gd, water_prop.kT, &potential, QuadraticLineMinimizer));
      retval += test_minimizer("PreconditionedConjugateGradient", pd, 100, 1e-5);

      double peak = peak_memory()/1024.0/1024;
      reset_peak_memory();
      printf("Peak memory use was %g\n", peak);
      if (peak > 7 || peak < 6) {
        printf("FAIL: Peak memory use was not as expected PreconditionedConjugateGradient...\n");
        retval++;
      }
    }

    //potential = external_potential + mu*VectorXd::Ones(gd.NxNyNz);
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    density.epsNativeSlice("cavity-density.eps", Cartesian(2*rmax,0,0), Cartesian(0,2*rmax,0), Cartesian(-rmax,-rmax,0));
    density.epsRadial1d("cavity-radial-density.eps", 0, rmax, nliquid, R, "Density scaled by nliquid");

    Grid grad(gd);
    grad.setZero();
    ff.integralgrad(water_prop.kT, potential, &grad);
 
    retval += constrain(constraint, f0wb).run_finite_difference_test("white bear functional",
                                                                     water_prop.kT, density, &grad);
    retval += constrain(constraint, f0rf).run_finite_difference_test("rosenfeld functional",
                                                                     water_prop.kT, density, &grad);
  }

  double peak = peak_memory()/1024.0/1024;
  printf("Peak memory use was %g\n", peak);
  if (peak > 17 || peak < 16) {
    printf("FAIL: Peak memory use was not as expected...\n");
    retval++;
  }

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
