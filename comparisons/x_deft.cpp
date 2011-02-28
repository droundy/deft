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
  double kT = water_prop.epsilon_dispersion*1.2;
  Functional f = Xassociation(R,
                              water_prop.epsilonAB, water_prop.kappaAB,
                              water_prop.epsilon_dispersion,
                              water_prop.lambda_dispersion);
  for (double eta=0.03125; eta<=0.5; eta += 0.03125) {
    double n = eta/(4*M_PI*R*R*R/3);
    // Energy units in the vrpack code are water_prop.epsilon_dispersion
    printf("%17.12f%17.12f\n", eta, f(kT, n));
  }
}
