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

//const double kT = water_prop.kT;
const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin
const double kT = kB*293; // Room temperature

int main(int, char **) {
  FILE *o = fopen("paper/figs/entropy.dat", "w");

  Functional f = SaftFluidSlow(water_prop.lengthscale, kT,
                               water_prop.epsilonAB, water_prop.kappaAB,
                               water_prop.epsilon_dispersion,
                               water_prop.lambda_dispersion, 0);
  double mu_satp = find_chemical_potential(f, kT,
                                           water_prop.liquid_density);
  f = SaftFluidSlow(water_prop.lengthscale, kT,
                    water_prop.epsilonAB, water_prop.kappaAB,
                    water_prop.epsilon_dispersion,
                    water_prop.lambda_dispersion, mu_satp);
  const double nl = saturated_liquid(f, kT);
  printf("Saturated liquid density turns out to be: %g\n", nl);
  const double nv = coexisting_vapor_density(f, kT, nl);
  printf("Vapor density turns out to be: %g\n", nv);
  
  Functional S = SaftEntropy(water_prop.lengthscale, kT, water_prop.epsilonAB, water_prop.kappaAB,
			       water_prop.epsilon_dispersion, water_prop.lambda_dispersion);
  //double mu = find_chemical_potential(f, kT, nl);
  for (double dens=0.01*nv; dens<=1.2*nl; dens *= 1.01) {
    double V = -kT*log(dens);
    //double Vl = -kT*log(nl);
    fprintf(o, "%g\t%g\t%g\t%g\n", dens, f(V), S(V), f(V) - kT*S(V)); //Prints n, F, S, U to data file
  }
  fclose(o);
}
