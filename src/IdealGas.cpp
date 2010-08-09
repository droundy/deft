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

#include "Functionals.h"
#include <stdio.h>
#include <math.h>

class IdealGasType : public FunctionalInterface {
public:
  IdealGasType(const GridDescription &g, double temp)
    : gd(g), T(temp) {}

  double energy(double data) const;
  double energy(const VectorXd &data) const;
  void grad(const VectorXd &data,
            VectorXd *, VectorXd *pgrad = 0) const;

  void  print_summary(const char *prefix, const VectorXd &data) const {
    printf("%sIdealGas energy = %g\n", prefix, (*this)(data));
  }
private:
  GridDescription gd;
  double T; // temperature
};

const double min_log_arg = 1e-10;
const double slope = log(min_log_arg);
const double min_e = min_log_arg*log(min_log_arg) - min_log_arg;

double IdealGasType::energy(const VectorXd &data) const {
  double e = 0;
  for (int i=0; i < gd.NxNyNz; i++) {
    const double n = data[i];
    if (n != n) {
      printf("data[%d] == %g\n", i, n);
      assert(n == n); // check for NaN
    }
    if (n > min_log_arg) {
      e += n*log(n) - n;
      if (isinf(n)) return INFINITY; // Infinite density gives infinite energy.
      assert(e == e); // check for NaN
    } else {
      e += (n-min_log_arg)*slope + min_e;
    }
  }
  return T*e*gd.Lat.volume()/gd.NxNyNz;
}

double IdealGasType::energy(double n) const {
  double e;
  if (n != n) {
    printf("n == %g\n", n);
    assert(n == n); // check for NaN
  }
  if (n > min_log_arg) {
    e = n*log(n) - n;
    if (isinf(n)) return INFINITY; // Infinite density gives infinite energy.
    assert(e == e); // check for NaN
  } else {
    e = (n-min_log_arg)*slope + min_e;
  }
  return T*e*gd.Lat.volume()/gd.NxNyNz;
}

void IdealGasType::grad(const VectorXd &n,
                    VectorXd *g_ptr, VectorXd *pg_ptr) const {
  VectorXd &g = *g_ptr;

  const double TdV = T*gd.Lat.volume()/gd.NxNyNz;
  for (int i=0; i < gd.NxNyNz; i++) {
    if (n[i] > min_log_arg) {
      g[i] += TdV*log(n[i]);
    } else {
      g[i] += TdV*slope;
    }
  }
  if (pg_ptr) *pg_ptr += g;
}

Functional IdealGas(const GridDescription &g, double temperature) {
  return Functional(new IdealGasType(g, temperature));
}
