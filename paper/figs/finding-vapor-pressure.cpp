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

#define NUMT 2

int main(int, char **) {
  //const double kT = water_prop.kT;
  const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin

  FILE *o = fopen("paper/figs/finding-vapor-pressure.dat", "w");

  double mu_satp;
  {
    Functional f = SaftFluid(water_prop.lengthscale,
                             water_prop.epsilonAB, water_prop.kappaAB,
                             water_prop.epsilon_dispersion,
                             water_prop.lambda_dispersion, water_prop.length_scaling, 0);
    mu_satp = find_chemical_potential(f, water_prop.kT,
                                      water_prop.liquid_density);
  }

  Functional f = SaftFluid(water_prop.lengthscale,
                           water_prop.epsilonAB, water_prop.kappaAB,
                           water_prop.epsilon_dispersion,
                           water_prop.lambda_dispersion, water_prop.length_scaling, mu_satp);

  //double Temperatures[NUMT] = { water_prop.kT/kB };
  double Temperatures[NUMT] = { water_prop.kT/kB, 693 };
  double mu[NUMT], nl[NUMT], nv[NUMT];

  double nvmin=1e100, nlmax=0;

  for (int t=0; t < NUMT; t++) {
    double kT = kB*Temperatures[t];

    mu[t] = 0; //find_chemical_potential(fs[t], kT, water_prop.vapor_density);
    saturated_liquid_vapor(f, kT, 1e-14, 0.0017, 0.0055, &nl[t], &nv[t], &mu[t], 1e-5);
    //printf("Saturated liquid density at %gK turns out to be: %g\n", Temperatures[t], nl[t]);
    //printf("Vapor density turns out to be: %g\n", nv[t]);

    //double mu2 = find_chemical_potential(f, kT, nl[t]);
    //printf("Chemical potential comparison: %g vs %g differ by %g\n", mu[t], mu2, mu[t] - mu2);

    if (nv[t] < nvmin) nvmin = nv[t];
    if (nl[t] > nlmax) nlmax = nl[t];
  }
  
  for (double dens=0.01*nvmin; dens<=1.2*nlmax; dens *= 1.01) {
    fprintf(o, "%g", dens);
    for (int t=0; t<NUMT; t++) {
      const double kT = kB*Temperatures[t];
      const double V = -kT*log(dens);
      const double Vl = -kT*log(nl[t]);
      fprintf(o, "\t%g\t%g", f(kT, V), f(kT, Vl) - (dens-nl[t])*mu[t]);
    }
    fprintf(o, "\n");
    fflush(o);
  }
  fclose(o);
}
