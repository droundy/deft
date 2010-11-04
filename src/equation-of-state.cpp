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

#include "utilities.h"
#include "Functional.h"
#include <stdio.h>

double find_density(Functional f, double nmin, double nmax) {
  double nbest = nmin;
  double ebest = f(nmin);
  const double factor = 1.1;
  for (double n = nmin; n<=nmax; n *= factor) {
    double en = f(n);
    if (en < ebest) {
      ebest = en;
      nbest = n;
    }
  }
  double nlo = nbest/factor, nhi = nbest*factor;
  while ((nhi - nlo)/nbest > 1e-7) {
    if (nbest < 0.5*(nhi+nlo)) {
      double ntry = 0.3*nlo + 0.7*nhi;
      double etry = f(ntry);
      if (etry < ebest) {
        nlo = nbest;
        nbest = ntry;
        ebest = etry;
      } else {
        nhi = ntry;
      }
    } else {
      double ntry = 0.7*nlo + 0.3*nhi;
      double etry = f(ntry);
      if (etry < ebest) {
        nhi = nbest;
        nbest = ntry;
        ebest = etry;
      } else {
        nlo = ntry;
      }
    }
  }
  return nbest;
}
