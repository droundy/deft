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

class Choose : public FieldFunctionalInterface {
public:
  Choose(double c, const FieldFunctional &fa, const FieldFunctional &fb) : cut(c), flow(fa), fhigh(fb) {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    VectorXd out(data);
    const int N = data.cols()*data.rows();
    for (int i=0; i<N; i++) {
      double x = data[i];
      out[i] = (x<cut) ? flow(x) : fhigh(x);
    }
    return out;
  }
  double transform(double x) const {
    return (x<cut) ? flow(x) : fhigh(x);
  }
  double grad(double x) const {
    return (x<cut) ? flow.grad(x) : fhigh.grad(x);
  }

  FieldFunctional grad(const FieldFunctional &ingrad) const {
    return ingrad*choose(cut, flow.grad(FieldFunctional(1)), fhigh.grad(FieldFunctional(1)));
  }
  void grad(const GridDescription &, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    const int N = data.cols()*data.rows();
    if (outpgrad) {
      for (int i=0; i<N; i++) {
        double x = data[i];
        if (x<cut) {
          (*outgrad)[i] += ingrad[i]*flow.grad(x);
          (*outpgrad)[i] += ingrad[i]*flow.grad(x);
        } else {
          (*outgrad)[i] += ingrad[i]*fhigh.grad(x);
           (*outpgrad)[i] += ingrad[i]*fhigh.grad(x);
        }
      }
    } else {
      for (int i=0; i<N; i++) {
        double x = data[i];
        (*outgrad)[i] += ingrad[i] * ((x<cut) ? flow.grad(x) : fhigh.grad(x));
      }
    }
  }
private:
  double cut;
  FieldFunctional flow, fhigh;
};

FieldFunctional choose(double cut, const FieldFunctional &lower, const FieldFunctional &higher) {
  return FieldFunctional(new Choose(cut, lower, higher));
}

static const double min_log_arg = 1e-90;
static const double slope = log(min_log_arg);
static const double min_e = min_log_arg*log(min_log_arg) - min_log_arg;

FieldFunctional IdealGas(double T) {
  FieldFunctional n = Identity().set_name("n");
  return choose(min_log_arg,  T*((n-min_log_arg)*slope + min_e), T*(n*log(n)-n)).set_name("ideal gas");
}
