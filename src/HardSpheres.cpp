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

/*
class HardSpheresType : public FunctionalInterface {
public:
  HardSpheresType(double rad, double temp) : R(rad), T(temp) {}

  double energy(double data) const;
  double energy(const GridDescription &gd, const VectorXd &data) const;
  void grad(const GridDescription &gd, const VectorXd &data,
            VectorXd *, VectorXd *pgrad = 0) const;

  void  print_summary(const char *prefix, double last_energy) const {
    printf("%sHardSpheres energy = %g\n", prefix, last_energy);
  }
private:
  double R, T; // temperature
};

const double min_log_arg = 1e-12;
const double slope = log(min_log_arg);
const double min_e = min_log_arg*log(min_log_arg) - min_log_arg;

double HardSpheresType::energy(const GridDescription &gd, const VectorXd &data) const {
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
  return T*e*gd.dvolume;
}

double HardSpheresType::energy(double n) const {
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
  return T*e;
}

void HardSpheresType::grad(const GridDescription &gd, const VectorXd &n,
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
*/

Functional HardSpheres(double radius, double temperature) {
  FieldFunctional n0 = Identity();
  FieldFunctional n3 = StepConvolve(radius);
  FieldFunctional one_minus_n3 = 1 - StepConvolve(radius);
  FieldFunctional n2 = ShellConvolve(radius);
  FieldFunctional n2x = xShellConvolve(radius);
  FieldFunctional n2y = yShellConvolve(radius);
  FieldFunctional n2z = zShellConvolve(radius);
  FieldFunctional nTxx = xxShellConvolve(radius);
  FieldFunctional nTyy = yyShellConvolve(radius);
  FieldFunctional nTzz = zzShellConvolve(radius);
  FieldFunctional nTxy = xyShellConvolve(radius);
  FieldFunctional nTyz = yzShellConvolve(radius);
  FieldFunctional nTzx = zxShellConvolve(radius*temperature);
  FieldFunctional phi1 = -1*n0*log(one_minus_n3);
  const double four_pi_r2 = (4*M_PI)*(radius*radius);
  // n1 is n2/(four_pi_r2)
  FieldFunctional phi2 = (n2*n2 - n2x*n2x - n2y*n2y - n2z*n2z)/(four_pi_r2*one_minus_n3);
  return temperature*integrate(phi1 + phi2);
    //integrate(phi2);
}
