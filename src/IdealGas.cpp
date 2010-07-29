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

#include "IdealGas.h"
#include <stdio.h>
#include <math.h>

void IdealGas::print_summary() const {
  printf("Ideal gas summary\n");
}

const double min_log_arg = 1e-10;
const double slope = log(min_log_arg);
const double min_e = min_log_arg*log(min_log_arg) - min_log_arg;

double IdealGas::operator()(const VectorXd &data) const {
  double e = 0;
  for (int i=0; i < gd.NxNyNz; i++) {
    const double n = data[i];
    assert(n == n); // check for NaN
    if (n > min_log_arg) {
      e += n*log(n) - n;
      assert(e == e); // check for NaN
    } else {
      e += (n-min_log_arg)*slope + min_e;
    }
  }
  return T*e;
}

void IdealGas::grad(const VectorXd &n,
                    VectorXd *g_ptr, VectorXd *pg_ptr) const {
  VectorXd &g = *g_ptr;

  for (int i=0; i < gd.NxNyNz; i++) {
    if (n[i] > min_log_arg) {
      g[i] += T*log(n[i]);
    } else {
      g[i] += T*slope;
    }
  }
  if (pg_ptr) *pg_ptr += g;
}
