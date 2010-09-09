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

static const double min_log_arg = 1e-90;
static const double slope = log(min_log_arg);
static const double min_e = min_log_arg*log(min_log_arg) - min_log_arg;

class IdealGasType : public FieldFunctionalInterface {
public:
  IdealGasType(double temp) : T(temp) {}


  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    const int N = data.cols()*data.rows();
    VectorXd out(data);
    for (int i=0; i<N; i++) {
      const double n = out[i];
      if (isnan(n)) {
        printf("data[%d] == %g\n", i, n);
        assert(!isnan(n));
      }
      if (n > min_log_arg) {
        if (isinf(n)) {
          out[i] = INFINITY; // Infinite density gives infinite energy.
        } else {
          out[i] = T*(n*log(n) - n);
          assert(!isnan(log(n)));
        }
      } else {
        out[i] = T*((n-min_log_arg)*slope + min_e);
      }
    }
    return out;
  }
  double transform(double n) const {
    assert(!isnan(n));
    if (n > min_log_arg) {
      if (isinf(n))
        return INFINITY; // Infinite density gives infinite energy.
      else
        return T*(n*log(n) - n);
    } else
      return T*((n-min_log_arg)*slope + min_e);
  }

  void grad(const GridDescription &, const VectorXd &n, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    const int N = n.cols()*n.rows();
    if (outpgrad)
      for (int i=0; i<N; i++) {
        if (n[i] > min_log_arg) {
          const double X = T*ingrad[i]*log(n[i]);
          (*outgrad)[i] += X;
          (*outpgrad)[i] += X;
        } else {
          const double X = T*ingrad[i]*slope;
          (*outgrad)[i] += X;
          (*outpgrad)[i] += X;
        }
      }
    else
      for (int i=0; i<N; i++) {
        if (n[i] > min_log_arg)
          (*outgrad)[i] += T*ingrad[i]*log(n[i]);
        else
          (*outgrad)[i] += T*ingrad[i]*slope;
      }
  }
private:
  double T; // temperature
};

FieldFunctional IdealGas(double temperature) {
  return FieldFunctional(new IdealGasType(temperature), "ideal gas energy");
}
