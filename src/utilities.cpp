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

#include "utilities.h"
#include "Grid.h"
#include "Functionals.h"
#include <stdio.h>

LiquidProperties water_prop = {
  3, // approx size of water in bohr
  4.9388942e-3, // density of liquid water
  1.14e-7, // vapor density of water
  1e-3 // room temperature in Hartree
};

double surface_tension(Minimizer min, FieldFunctional f, LiquidProperties prop, bool verbose) {
  int numptspersize = 100;
  int size = 16;
  Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,size*prop.lengthscale));
  GridDescription gd(lat, 1, 1, numptspersize*size);
  Grid potential(gd);

  // Set the density to range from vapor to liquid
  const double Veff_liquid = -prop.kT*log(prop.liquid_density);
  const double Veff_gas = -prop.kT*log(prop.vapor_density);
  for (int i=0; i<gd.NxNyNz/2; i++) potential[i] = Veff_gas;
  for (int i=gd.NxNyNz/2; i<gd.NxNyNz; i++) potential[i] = Veff_liquid;

  FieldFunctional n = EffectivePotentialToDensity(prop.kT);
  FieldFunctional f0 = f(n);
  f0.run_finite_difference_test("f0", potential);
  min.minimize(f0, gd, &potential);
  while (min.improve_energy(verbose))
    if (verbose) {
      printf("Working on liberated interface...\n");
      fflush(stdout);
    }
  const double Einterface = f0.integral(potential);
  double Ninterface = 0;
  {
    Grid density(gd, EffectivePotentialToDensity(prop.kT)(gd, potential));
    for (int i=0;i<gd.NxNyNz;i++) Ninterface += density[i]*gd.dvolume;
  }
  if (verbose) printf("Got interface energy of %g.\n", Einterface);
  
  for (int i=0; i<gd.NxNyNz; i++) potential[i] = Veff_gas;
  min.minimize(f0, gd, &potential);
  while (min.improve_energy(verbose))
    if (verbose) {
      printf("Working on gas...\n");
      fflush(stdout);
    }
  const double Egas = f0.integral(potential);
  double Ngas = 0;
  {
    Grid density(gd, EffectivePotentialToDensity(prop.kT)(gd, potential));
    for (int i=0;i<gd.NxNyNz;i++) Ngas += density[i]*gd.dvolume;
  }
  
  for (int i=0; i<gd.NxNyNz; i++) potential[i] = Veff_liquid;
  min.minimize(f0, gd, &potential);
  while (min.improve_energy(verbose))
    if (verbose) {
      printf("Working on liquid...\n");
      fflush(stdout);
    }
  const double Eliquid = f0.integral(potential);
  double Nliquid = 0;
  {
    Grid density(gd, EffectivePotentialToDensity(prop.kT)(gd, potential));
    for (int i=0;i<gd.NxNyNz;i++) Nliquid += density[i]*gd.dvolume;
  }
  
  const double X = Ninterface/Nliquid; // Fraction of volume effectively filled with liquid.
  const double surface_tension = (Einterface - Eliquid*X - Egas*(1-X))/2;
  if (verbose) {
    printf("\n\n");
    printf("interface energy is %.15g\n", Einterface);
    printf("gas energy is %.15g\n", Egas);
    printf("liquid energy is %.15g\n", Eliquid);
    printf("Ninterface/liquid/gas = %g/%g/%g\n", Ninterface, Nliquid, Ngas);
    printf("X is %g\n", X);
    printf("surface tension is %.10g\n", surface_tension);
  }
  return surface_tension;
}

long Eigen::djr_memused = 0;
long Eigen::djr_mempeak = 0;

long peak_memory() {
  return Eigen::djr_mempeak;
}

void reset_peak_memory() {
  Eigen::djr_mempeak = Eigen::djr_memused;
}
