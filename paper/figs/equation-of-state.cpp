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
#include "equation-of-state.h"
#include "steam-table.h"

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  printf("%s took %g seconds\n", name, (t-last_time)/double(CLOCKS_PER_SEC));
  last_time = t;
}

const double kB = 3.16681539628059e-6; // Boltzmann's constant in Hartree/Kelvin

int main(int, char **) {
  FILE *o = fopen("paper/figs/equation-of-state.dat", "w");
  FILE *experiment = fopen("paper/figs/experimental-equation-of-state.dat", "w");
  int imax=0;
  while (temperatures_kelvin[imax]) imax++;
  took("Counting the temperatures");
  Functional f = SaftFluidSlow(water_prop.lengthscale,
                               water_prop.epsilonAB, water_prop.kappaAB,
                               water_prop.epsilon_dispersion,
                               water_prop.lambda_dispersion, 0);
  for (int i=0; i<imax; i+=3) {
    printf("Working on equation of state at %g Kelvin...\n", temperatures_kelvin[i]);
    double kT = kB*temperatures_kelvin[i];
    double nl = saturated_liquid(f, kT, 0.0015, 0.0055);
    took("Finding saturated liquid density");
    double nv = coexisting_vapor_density(f, kT, nl);
    took("Finding coesisting vapor density");
    double pv = pressure(f, kT, nv);
    took("Finding pressure");
      
    fprintf(o, "%g\t%g\t%g\t%g\n",
            temperatures_kelvin[i], pv, nl, nv);
    fflush(o); // FOR DEBUGGING
    fprintf(experiment, "%g\t%g\t%g\t%g\t%g\n",
            temperatures_kelvin[i], water_vapor_pressure[i],
            water_saturation_liquid[i], water_vapor_density[i],
            water_saturated_surface_tension[i]);
  }

  for (double T=660; T<700; T += 5) {
    printf("Working on bonus equation of state at %g Kelvin...\n", T);
    double kT = kB*T;
    const double vapor_liquid_ratio = 0.95;
    double nl = saturated_liquid(f, kT, 0.001, 0.005, vapor_liquid_ratio);
    took("Finding saturated liquid density");
    double nv = coexisting_vapor_density(f, kT, nl, vapor_liquid_ratio);
    took("Finding coesisting vapor density");
    double pv = pressure(f, kT, nv);
    took("Finding pressure");
      
    fprintf(o, "%g\t%g\t%g\t%g\n", T, pv, nl, nv);
    fflush(o); // FOR DEBUGGING
  }

  fclose(o);
  fclose(experiment);
}
