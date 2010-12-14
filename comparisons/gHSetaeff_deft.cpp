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
  Functional gHS = gHScarnahan(Identity(), R);
  Functional eta_eff = eta_effective(Identity(), water_prop.lambda_dispersion);
  for (double eta=0.0; eta<=0.5; eta += 0.03125) {
    printf("%15.10f%15.10f\n", eta, gHS(eta_eff(eta)));
  }
}
