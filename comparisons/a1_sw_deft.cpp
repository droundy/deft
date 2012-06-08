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
#include "Functionals.h"
#include "Testables.h"
#include "equation-of-state.h"

int main(int, char **) {
  double R = water_prop.lengthscale;
  Functional a1 = DispersionSAFTa1(R, water_prop.epsilon_dispersion,
                                   water_prop.lambda_dispersion, water_prop.length_scaling);
  for (double eta=0.0; eta<=0.5; eta += 0.03125) {
    double n = eta/(4*M_PI*R*R*R/3);
    // Energy units in the vrpack code are water_prop.epsilon_dispersion?
    double a1n = a1(0, n)/water_prop.epsilon_dispersion;
    // The following is a hokey trick to deal with optimizations like
    // -ffast-math that may turn -0.0 into 0.0.
    if (fabs(a1n) < 1e-12) printf("%17.12f  -0.000000000000\n", eta);
    else printf("%17.12f%17.12f\n", eta, a1n);
  }
}
