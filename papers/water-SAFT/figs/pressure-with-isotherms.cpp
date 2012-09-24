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
#include "equation-of-state.h"

//const double kT = water_prop.kT;
const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin
//const double kT = kB*293; // Temperature

int main(int, char **) {
  
  FILE *o = fopen("papers/water-SAFT/figs/pressure-with-isotherms.dat", "w");
  
  Functional f = OfEffectivePotential(SaftFluid2(water_prop.lengthscale,
                                                water_prop.epsilonAB, water_prop.kappaAB,
                                                water_prop.epsilon_dispersion,
                                                water_prop.lambda_dispersion, water_prop.length_scaling, 0));
  for (double dens=0.00001; dens<=0.0055; dens *= 1.002) {
    fprintf(o, "%g", dens);

    for (double kT=kB*298; kT<=kB*648; kT+=50*kB) {
      double mu = 0, nl = 0, nv = 0;
      saturated_liquid_vapor(f, kT, 1e-14, 0.0017, 0.0055, &nl, &nv, &mu, 1e-6);
      double p = 0;
      if (dens <= nv || dens >= nl) p = pressure(f, kT, dens);
      else p = pressure(f, kT, nv);
      //printf("Pressure = %g\n", p); //DEBUGGING
      fprintf(o, "\t%g\t%g", kT, p); //Prints kT, pressure, to data file
      //fflush(o);
  }
    fprintf(o, "\n");
  }
    fclose(o);
}
