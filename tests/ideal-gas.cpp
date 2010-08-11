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

int main() {
  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  int resolution = 20;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid density(gd);
  const double kT = 1e-3; // room temperature in Hartree
  density = 10*(-500*density.r2()).cwise().exp()
    + 1e-7*VectorXd::Ones(gd.NxNyNz);
  Functional ig = IdealGas(gd, kT);
  if (ig.run_finite_difference_test("ideal gas", density)) {
    printf("Finite difference test failed!!!\n");
    return 1;
  } else {
    printf("Finite difference test passes.\n");
  }

  printf("\nNow let's try this with an effective potential...\n");
  Grid potential(gd);
  potential = 1e-2*((-50*density.r2()).cwise().exp())
    + -2e-3*VectorXd::Ones(gd.NxNyNz);
  //potential.epsNativeSlice("potential.eps", Cartesian(1,0,0),
  //                         Cartesian(0,1,0), Cartesian(0,0,0));
  //density = EffectivePotentialToDensity(kT)(potential);
  //density.epsNativeSlice("dens.eps", Cartesian(1,0,0),
  //                         Cartesian(0,1,0), Cartesian(0,0,0));
  Functional ig2 = compose(IdealGas(gd,kT), EffectivePotentialToDensity(kT));
  density *= 1e-5;
  if (ig2.run_finite_difference_test("ideal gas", potential)) {
    printf("Finite difference test failed!!!\n");
    return 1;
  } else {
    printf("Finite difference test passes.\n");
  }
  return 0;
}
