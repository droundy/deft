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
#include "ExternalPotential.h"
#include "Downhill.h"

int main() {
  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  int resolution = 5;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid external_potential(gd);
  external_potential = 1e-1*external_potential.r2(); // a harmonic trap...
  Grid potential(gd);
  const double kT = 1e-3; // room temperature in Hartree
  const double ngas = 1e-5; // vapor density of water
  const double mu = -kT*log(ngas);
  potential = +1e-4*((-10*potential.r2()).cwise().exp())
    + 1.04*mu*VectorXd::Ones(gd.NxNyNz);
  //potential.epsNativeSlice("potential.eps", Cartesian(1,0,0),
  //                         Cartesian(0,1,0), Cartesian(0,0,0));
  counted_ptr<FunctionalSum> ig_and_mu =
    counted_ptr<Functional>(new IdealGas(gd,kT)) +
    counted_ptr<Functional>(new ChemicalPotential(gd, mu)) +
    counted_ptr<Functional>(new ExternalPotential(external_potential));
  counted_ptr<Functional> f =
    compose(counted_ptr<Functional>(ig_and_mu),
            counted_ptr<FieldFunctional>(new EffectivePotentialToDensity(kT)));
  Grid old_potential(potential);
  Downhill min(f, &potential);
  for (int i=0;i<30000;i++) { // using Downhill is *extremely* painful!
    min.improve_energy();
    min.print_info(i);
  }
  Grid density(gd);
  density = EffectivePotentialToDensity(kT)(potential);
  Grid expected_density(gd);
  expected_density = EffectivePotentialToDensity(kT)(external_potential + mu*VectorXd::Ones(gd.NxNyNz));
  double err2 = 0;
  for (int i=0;i<gd.NxNyNz;i++) {
    err2 += (density[i]-expected_density[i])*(density[i]-expected_density[i]);
    //err2 += (potential[i] - external_potential[i] - mu)*(potential[i] - external_potential[i] - mu);
  }
  err2 /= gd.NxNyNz;
  printf("rms error = %g\n", sqrt(err2));
  for (int i=0;i<gd.NxNyNz;i++) {
    //if (fabs(potential[i] - external_potential[i] - mu) > fabs(1e-7*(external_potential[i]+mu))) {
    if (fabs(density[i] - expected_density[i]) > fabs(1e-4*ngas)) {
      printf("Oh no, the error in potential is %g out of %g at %d!\n",
             potential[i]-external_potential[i]-mu, external_potential[i]+mu, i);
      printf("Oh no, the error in density is %g out of %g at %d!\n",
             density[i] - expected_density[i], expected_density[i], i);
      return 1;
    }
  }
}
