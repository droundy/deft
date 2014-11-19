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
#include "utilities.h"
#include "handymath.h"

Functional SoftFluid(double sigma, double epsilon, double mu);
Functional HardFluid(double radius, double mu);

int main(int, char **) {
  double radius = 1.0;
  double sigma = radius*pow(2,5.0/6.0);
  FILE *out = fopen("papers/fuzzy-fmt/figs/homogeneous.dat", "w");
  const double temps[] = { 0.0, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0};
  fprintf(out, "# eta");
  for (unsigned int i = 0; i<sizeof(temps)/sizeof(temps[0]); i++) {
    fprintf(out, "\tp(kT=%g)/nkT", temps[i]);
  }
  fprintf(out, "\n0");
  for (unsigned int i = 0; i<sizeof(temps)/sizeof(temps[0]); i++) {
    fprintf(out, "\t%g", temps[i]);
  }
  fprintf(out, "\n");
  for (double eta = 0.001; eta <= 0.7; eta *= 1.001) {
    fprintf(out, "%g", eta);
    for (unsigned int i = 0; i<sizeof(temps)/sizeof(temps[0]); i++) {
      const double temp = temps[i];
      Functional f = HardFluid(radius,0);
      if (temp > 0) f = SoftFluid(sigma, 1, 0);
      double usekT = temp;
      if (temp == 0) usekT = 1.0;
      const double n = eta/(4*M_PI/3);
      fprintf(out, "\t%g", pressure(OfEffectivePotential(f), usekT, n)/(n*usekT));
    }
    fprintf(out, "\n");
  }
  fclose(out);
  return 0;
}
