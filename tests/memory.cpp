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
#include "equation-of-state.h"
#include "utilities.h"

int retval = 0;
clock_t last_time = clock();
double reference_time = 1;

double check_peak(const char *name, const char *name2, double peakmin, double peakmax) {
  printf("\n===> Testing %s of %s <===\n", name, name2);
  
  double cputime = (clock()-last_time)/double(CLOCKS_PER_SEC);
  double peak = peak_memory()/1024.0/1024;
  printf("CPU time is %g s (or %gx)\n", cputime, cputime/reference_time);
  printf("Peak memory use is %g M\n", peak);
  if (peak < peakmin) {
    printf("FAIL: Peak memory use of %s %s should be at least %g!\n", name, name2, peakmin);
    retval++;
  }
  if (peak > peakmax) {
    printf("FAIL: Peak memory use of %s %s should be under %g!\n", name, name2, peakmax);
    retval++;
  }
  reset_peak_memory();
  last_time = clock();
  return cputime;
}

const double R = 2.7;
const double rcav = R+R; // 11.8*R+R;

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

void check_a_functional(const char *name, Functional f, const Grid &x,
                        double memE, double memG, double memP, double memPonly) {
  const double kT = water_prop.kT; // room temperature in Hartree

  printf("\n************");
  for (unsigned i=0;i<strlen(name) + 4;i++) printf("*");
  printf("\n* Working on %s *\n", name);
  for (unsigned i=0;i<strlen(name) + 4;i++) printf("*");
  printf("************\n\n");
  fflush(stdout);
  reset_peak_memory();
  last_time = clock();

  f.integral(kT, x);
  //printf("\n\nEnergy of %s is %g\n", name, f.integral(x));

  check_peak("Energy", name, memE, memE + 1);

  Grid mygrad(x);
  mygrad.setZero();
  f.integralgrad(kT, x, &mygrad);
  //printf("Grad of %s is: %g\n", name, mygrad.norm());

  check_peak("Gradient", name, memG, memG+1);
  
  {
    Grid mypgrad(x);
    mygrad.setZero();
    mypgrad.setZero();
    f.integralgrad(kT, x, &mygrad, &mypgrad);
  }
  check_peak("Gradient and preconditioned gradient", name, memP, memP+1);

  f.integralpgrad(kT, x, &mygrad);
  check_peak("Preconditioned gradient", name, memPonly, memPonly+1);

}

int main(int, char **argv) {
  const double kT = water_prop.kT; // room temperature in Hartree
  const double eta_one = 3.0/(4*M_PI*R*R*R);
  const double nliquid = 0.324*eta_one;
  Functional n = EffectivePotentialToDensity();
  const double mu = find_chemical_potential(HardSpheres(R)(n) + IdealGasOfVeff, kT, nliquid);

  // Here we set up the lattice.
  const double rmax = rcav*2;
  Lattice lat(Cartesian(0,rmax,rmax), Cartesian(rmax,0,rmax), Cartesian(rmax,rmax,0));
  //Lattice lat(Cartesian(1.4*rmax,0,0), Cartesian(0,1.4*rmax,0), Cartesian(0,0,1.4*rmax));
  GridDescription gd(lat, 0.2);

  last_time = clock();
  Grid external_potential(gd);
  // Do some pointless stuff so we can get some sort of gauge as to
  // how fast this CPU is, for comparison with other tests.
  for (int i=0; i<10; i++) {
    // Do this more times, to get a more consistent result...
    external_potential = external_potential.Ones(gd.NxNyNz);
    external_potential = external_potential.cwise().exp();
    external_potential = 13*external_potential + 3*external_potential.cwise().square();
    external_potential.fft(); // compute and toss the fft...
  }
  // And now let's set the external_potential up as we'd like it.
  external_potential.Set(incavity);
  external_potential *= 1e9;
  reference_time = check_peak("Setting", "external potential", 7, 8);

  Grid constraint(gd);
  constraint.Set(notincavity);
  //Functional f1 = f0 + ExternalPotential(external_potential);
  Functional ff = constrain(constraint, IdealGasOfVeff + (HardSpheres(R) + ChemicalPotential(mu))(n));
  
  Grid potential(gd, external_potential + 0.005*VectorXd::Ones(gd.NxNyNz));

  check_a_functional("HardSpheres", ff, potential, 83, 104, 111, 108);

  ff = constrain(constraint, (HardSpheresFast(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff);
  check_a_functional("HardSphereFast", ff, potential, 69, 94, 101, 97);

  ff = constrain(constraint, (HardSpheresRF(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff);
  check_a_functional("HardSphereRF", ff, potential, 55, 87, 94, 90);

  ff = constrain(constraint, (HardSpheresRFFast(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff);
  check_a_functional("HardSphereRFFast", ff, potential, 48, 73, 80, 76);

  ff = constrain(constraint, (HardSpheresTarazona(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff);
  check_a_functional("HardSphereTarazona", ff, potential, 87, 108, 114, 111);

  ff = constrain(constraint, (HardSpheresTarazonaFast(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff);
  check_a_functional("HardSphereTarazonaFast", ff, potential, 69, 94, 101, 97);

  ff = constrain(constraint, (HardSpheresWBnotensor(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff);
  check_a_functional("HardSpheresWBnotensor", ff, potential, 55, 90, 97, 94);

  ff = constrain(constraint, (HardSpheresNoTensor(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff);
  check_a_functional("HardSphereNoTensor", ff, potential, 48, 73, 80, 76);

  ff = constrain(constraint, HardSphereGas(R, mu));
  check_a_functional("HardSphereGas", ff, potential, 66, 87, 94, 87);

  double eps = water_prop.epsilonAB;
  double kappa = water_prop.kappaAB;
  ff = SaftFluid(R, eps, kappa, water_prop.epsilon_dispersion, water_prop.lambda_dispersion, mu);
  check_a_functional("SaftFluid", ff, potential, 48, 87, 90, 87);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
