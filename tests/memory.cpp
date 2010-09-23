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

double check_peak(const char *name, double peakmin, double peakmax, double time_limit) {
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");
  
  double cputime = (clock()-last_time)/double(CLOCKS_PER_SEC);
  double peak = peak_memory()/1024.0/1024;
  printf("Peak memory use is %g M\n", peak);
  printf("CPU time is %g s (or %gx)\n", cputime, cputime/reference_time);
  if (peak < peakmin) {
    printf("FAIL: Peak memory use should be at least %g!\n", peakmin);
    retval++;
  }
  if (peak > peakmax) {
    printf("FAIL: Peak memory use should be under %g!\n", peakmax);
    retval++;
  }
  if (time_limit > 0 && cputime/reference_time > time_limit) {
    printf("FAIL: CPU time used should be under %gx!\n", time_limit);
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
  external_potential.Set(incavity);
  external_potential = external_potential.cwise().exp();
  external_potential = 13*external_potential.cwise().log() + 3*external_potential.cwise().square();
  external_potential.Set(incavity);
  external_potential *= 1e9;
  reference_time = check_peak("Setting external potential.", 7, 8, 0);

  Grid constraint(gd);
  constraint.Set(notincavity);
  //Functional f1 = f0 + ExternalPotential(external_potential);
  Functional n = EffectivePotentialToDensity(kT);
  Functional ff = constrain(constraint, (HardSpheres(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  
  Grid potential(gd, external_potential + 0.005*VectorXd::Ones(gd.NxNyNz));
  printf("Energy is %g\n", ff.integral(potential));

  check_peak("Energy of 3D cavity", 83, 84, 40);

  Grid mygrad(gd);
  mygrad.setZero();
  ff.integralgrad(potential, &mygrad);
  printf("Grad is: %g\n", mygrad.norm());

  check_peak("Gradient of 3D cavity", 100, 105, 250);

  ff = constrain(constraint, (HardSpheresFast(R, kT) + IdealGas(kT) + ChemicalPotential(mu))(n));
  printf("Energy new way is %g\n", ff.integral(potential));

  check_peak("Energy of 3D cavity (fast)", 69, 70, 10);

  mygrad.setZero();
  ff.integralgrad(potential, &mygrad);
  printf("Grad optimized is: %g\n", mygrad.norm());

  check_peak("Gradient of 3D cavity (fast)", 363, 364, 60);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
