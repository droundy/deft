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
#include "equation-of-state.h"

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  printf("%s took %g seconds\n", name, (t-last_time)/double(CLOCKS_PER_SEC));
  last_time = t;
}


int main(int, char **) {
  FILE *o = fopen("paper/figs/surface-tension.dat", "w");
  for (double T=273; T<=693; T+=25) {
    printf("Working on surface tension at %g Kelvin...\n", T);
    const double kB = 3.16681539628059e-6; // Boltzmann's constant in Hartree/Kelvin
    LiquidProperties prop = water_prop;
    prop.kT = kB*T;
    Functional f = SaftFluid(prop.lengthscale,
                             prop.epsilonAB, prop.kappaAB,
                             prop.epsilon_dispersion,
                             prop.lambda_dispersion, 0);
    saturated_liquid_properties(f, &prop);
    took("Finding bulk densities");
    double mu = find_chemical_potential(f, prop.kT, prop.liquid_density);
    f = SaftFluid(prop.lengthscale,
                  prop.epsilonAB, prop.kappaAB,
                  prop.epsilon_dispersion,
                  prop.lambda_dispersion, mu);
    char *plotname = (char *)malloc(1024);
    sprintf(plotname, "paper/figs/surface-%03g.dat", T);
    // Here we set up an unused lattice.
    Lattice lat(Cartesian(0.2,0,0), Cartesian(0,0.2,0), Cartesian(0,0,20));
    GridDescription gd(lat, 1, 1, 200);
    Grid foo(gd);
    Minimizer min = Precision(1e-7, ConjugateGradient(f, gd, prop.kT, &foo, QuadraticLineMinimizer));
    double st = surface_tension(min, f, prop, true, plotname);
    free(plotname);
    fprintf(o, "%g\t%g\n", T, st);
    if (T > 373 && T < 600) T += 25; // Use higher interval at higher temperatures
  }
  fclose(o);
}
