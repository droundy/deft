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
#include "ChemicalPotential.h"
#include "EffectivePotentialToDensity.h"
#include "Downhill.h"

int main() {
  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  int resolution = 5;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid potential(gd);
  const double kT = 1e-3; // room temperature in Hartree
  const double ngas = 1e-5; // vapor density of water
  const double mu = -kT*log(ngas);
  potential = +1e-4*((-10*potential.r2()).cwise().exp())
    + 1.04*mu*VectorXd::Ones(gd.NxNyNz);
  //potential.epsNativeSlice("potential.eps", Cartesian(1,0,0),
  //                         Cartesian(0,1,0), Cartesian(0,0,0));
  //density.epsNativeSlice("dens.eps", Cartesian(1,0,0),
  //                         Cartesian(0,1,0), Cartesian(0,0,0));
  Functional ig_and_mu = IdealGas(gd,kT) + ChemicalPotential(gd, mu);
  Functional f = compose(ig_and_mu, EffectivePotentialToDensity(kT));
  Grid old_potential(potential);
  Downhill min(f, &potential);
  for (int i=0;i<500;i++) {
    min.improve_energy();
    min.print_info(i);
  }
  double err2 = 0;
  for (int i=0;i<gd.NxNyNz;i++) {
    err2 += (potential[i] - mu)*(potential[i] - mu);
  }
  err2 /= gd.NxNyNz;
  printf("rms error = %g\n", sqrt(err2));
  for (int i=0;i<gd.NxNyNz;i++) {
    if (fabs(potential[i] - mu) > fabs(1e-10*mu)) {
      printf("Oh no, the error is %g out of %g at %d!\n", potential[i]-mu, mu, i);
      return 1;
    }
  }
}
