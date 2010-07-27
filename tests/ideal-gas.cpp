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
#include "IdealGas.h"
#include "EffectivePotentialToDensity.h"

int main() {
  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  int resolution = 20;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid density(gd);
  density = 10*(-500*density.r2()).cwise().exp()
    + 1e-7*VectorXd::Ones(gd.NxNyNz);
  IdealGas ig(gd, 300);
  if (ig.run_finite_difference_test("ideal gas", density)) {
    printf("Finite difference test passes.\n");
  } else {
    printf("Finite difference test failed!!!\n");
    return 1;
  }

  printf("\nNow let's verify it works with very low density!\n");
  density *= 1e-5;
  if (ig.run_finite_difference_test("ideal gas", density)) {
    printf("Finite difference test passes.\n");
  } else {
    printf("Finite difference test failed!!!\n");
    return 1;
  }

  printf("\nNow let's try this with an effective potential...\n");
  FunctionalComposition ig2 =
    *compose(counted_ptr<Functional>(new IdealGas(gd,300)),
            counted_ptr<FieldFunctional>(new EffectivePotentialToDensity(300)));
  density *= 1e-5;
  if (ig.run_finite_difference_test("ideal gas", density)) {
    printf("Finite difference test passes.\n");
  } else {
    printf("Finite difference test failed!!!\n");
    return 1;
  }
  return 0;
}
