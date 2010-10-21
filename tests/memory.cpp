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

int retval = 0;
clock_t last_time = clock();
double reference_time = 1;

double check_peak(const char *name, const char *name2, double peakmin, double peakmax, double time_limit, double fraccuracy=0.2) {
  printf("\n===> Testing %s of %s <===\n", name, name2);
  
  double cputime = (clock()-last_time)/double(CLOCKS_PER_SEC);
  double peak = peak_memory()/1024.0/1024;
  printf("CPU time is %g s (or %gx)\n", cputime, cputime/reference_time);
  if (time_limit > 0 && cputime/reference_time > time_limit*(1+fraccuracy)) {
    printf("FAIL: CPU time used by %s %s should be under %gx!\n", name, name2, time_limit*(1+fraccuracy));
    retval++;
  }
  if (cputime/reference_time < time_limit*(1-fraccuracy)) {
    printf("FAIL: CPU time used by %s %s should be at least %gx!\n", name, name2, time_limit*(1-fraccuracy));
    retval++;
  }
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
                        double memE, double timE,
                        double memG, double timG,
                        double memP, double timP,
                        double memPonly, double timPonly) {

  printf("\n************");
  for (unsigned i=0;i<strlen(name) + 4;i++) printf("*");
  printf("\n* Working on %s *\n", name);
  for (unsigned i=0;i<strlen(name) + 4;i++) printf("*");
  printf("************\n\n");
  fflush(stdout);
  reset_peak_memory();
  last_time = clock();

  f.integral(x);
  //printf("\n\nEnergy of %s is %g\n", name, f.integral(x));

  check_peak("Energy", name, memE, memE + 1, timE);

  Grid mygrad(x);
  mygrad.setZero();
  f.integralgrad(x, &mygrad);
  //printf("Grad of %s is: %g\n", name, mygrad.norm());

  check_peak("Gradient", name, memG, memG+1, timG);
  
  {
    Grid mypgrad(x);
    mygrad.setZero();
    mypgrad.setZero();
    f.integralgrad(x, &mygrad, &mypgrad);
  }
  check_peak("Gradient and preconditioned gradient", name, memP, memP+1, timP);

  f.integralpgrad(x, &mygrad);
  check_peak("Preconditioned gradient", name, memPonly, memPonly+1, timPonly);

}

int main(int, char **argv) {
  const double kT = water_prop.kT; // room temperature in Hartree
  const double eta_one = 3.0/(4*M_PI*R*R*R);
  const double nliquid = 0.324*eta_one;
  const double mu = -(HardSpheres(R, kT) + IdealGas(kT)).derive(nliquid);

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
  reference_time = check_peak("Setting", "external potential", 7, 8, 0);

  Grid constraint(gd);
  constraint.Set(notincavity);
  //Functional f1 = f0 + ExternalPotential(external_potential);
  Functional n = EffectivePotentialToDensity(kT);
  Functional ff = constrain(constraint, (HardSpheres(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  
  Grid potential(gd, external_potential + 0.005*VectorXd::Ones(gd.NxNyNz));

  check_a_functional("HardSpheres", ff, potential, 80, 10.5, 101, 70, 108, 70, 108, 70);

  ff = constrain(constraint, (HardSpheresFast(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  check_a_functional("HardSphereFast", ff, potential, 62, 2.0, 104, 15, 111, 15, 111, 15);

  ff = constrain(constraint, (HardSpheresRF(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  check_a_functional("HardSphereRF", ff, potential, 52, 3.5, 83, 17, 90, 17, 90, 17);

  ff = constrain(constraint, (HardSpheresRFFast(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  check_a_functional("HardSphereRFFast", ff, potential, 41, 1.2, 66, 3.8, 73, 3.8, 73, 3.8);

  ff = constrain(constraint, (HardSpheresTarazona(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  check_a_functional("HardSphereTarazona", ff, potential, 83, 10.0, 104, 62.1, 111, 62, 111, 62);

  ff = constrain(constraint, (HardSpheresTarazonaFast(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  check_a_functional("HardSphereTarazonaFast", ff, potential, 62, 2.0, 94, 10.2, 101, 10, 101, 10);

  ff = constrain(constraint, (HardSpheresWBnotensor(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  check_a_functional("HardSpheresWBnotensor", ff, potential, 52, 4.2, 87, 22, 94, 22, 94, 22);

  ff = constrain(constraint, (HardSpheresNoTensor(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  //check_a_functional("HardSphereNoTensor", ff, potential, 41, 1.2, 66, 5.0);
  check_a_functional("HardSphereNoTensor", ff, potential, 41, 1.2, 80, 5.0, 87, 5, 87, 5);

  ff = constrain(constraint, HardSphereGas(R, kT, mu));
  check_a_functional("HardSphereGas", ff, potential, 62, 2.0, 97, 13.8, 108, 52.6, 108, 52.6);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
