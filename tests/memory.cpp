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
const double eta_one = 3.0/(4*M_PI*R*R*R);
const double diameter_cubed = 1/(8*R*R*R);
const double nliquid = 0.324*eta_one;
const double mu = -integrate(HardSpheres(R, kT) + IdealGas(kT)).grad(nliquid);

// Here we set up the lattice.
const double rcav = R+R; // 11.8*R+R;
const double rmax = rcav*2;
Lattice lat(Cartesian(0,rmax,rmax), Cartesian(rmax,0,rmax), Cartesian(rmax,rmax,0));
//Lattice lat(Cartesian(1.4*rmax,0,0), Cartesian(0,1.4*rmax,0), Cartesian(0,0,1.4*rmax));
GridDescription gd(lat, 0.2);

// And the functional...
Functional f0 = integrate(HardSpheres(R, kT) + IdealGas(kT) + ChemicalPotential(mu));
Functional f0wb = integrate(HardSpheresWB(R, kT));
Functional f0rf = integrate(HardSpheresRF(R, kT));
FieldFunctional n = EffectivePotentialToDensity(kT);
Functional f = f0(n);

Grid external_potential(gd);
Grid potential(gd);
Functional ff;

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
  external_potential.Set(incavity);
  external_potential *= 1e9;
  external_potential.epsNativeSlice("external.eps", Cartesian(2*rmax,0,0), Cartesian(0,2*rmax,0), Cartesian(-rmax,-rmax,0));
  external_potential.epsRadial1d("external-radial.eps", 0, rmax, 1, R, "Good fun!");
  Grid constraint(gd);
  constraint.Set(notincavity);
  //Functional f1 = f0 + ExternalPotential(external_potential);
  Functional f1 = constrain(constraint, f0);
  ff = f1(n);
  
  int retval = 0;

  potential = external_potential + 0.005*VectorXd::Ones(gd.NxNyNz);
  printf("Energy is %g\n", ff(potential));

  double peak = peak_memory()/1024.0/1024;
  const double expected_peak = 70;
  printf("Peak memory use is %g M\n", peak);
  if (peak > expected_peak) {
    printf("FAIL: Peak memory use should be under %g!\n", expected_peak);
    retval++;
  }


  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
