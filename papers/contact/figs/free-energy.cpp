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
#include "ContactDensity.h"


int main(int, char **) { 
  FILE *o = fopen("papers/contact/figs/free-energy.dat", "w");
  assert(o);

  Functional f = HardSpheresWBnotensor(1.0);
  Functional eta = Identity();
  Functional fcar = (4*eta - 3*sqr(eta))/sqr(1-eta);
  double mykT = 1.0; // so we needn't divide by kT

  for (double eta=0.0001; eta<=0.3; eta *= 1.01) {
    double fc = fcar(mykT, eta);
    double n = eta/(4*M_PI/3);
    double nice = f(mykT, n)/n;
    fprintf(o, "%g\t%g\t%g\n", eta, fc, nice);
    fflush(o);
  }
  fclose(o);
}
