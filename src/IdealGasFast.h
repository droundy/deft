// -*- mode: C++; -*-

#pragma once

#include "Functional.h"

class IdealGasFast_type : public FunctionalInterface {
public:
  IdealGasFast_type(double kT_arg) : kT(kT_arg) { have_analytic_grad = false; }
  double transform(double) const {
    return 0;
  }
  double derive(double) const {
    return 0;
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &x) const {
    assert(&gd); // to avoid an unused parameter error
    assert(&x); // to avoid an unused parameter error
    return kT*choose(1e-90, -207.2326583694641*(x - 1e-90*VectorXd::Ones(gd.NxNyNz)) - 2.082326583694641e-88*VectorXd::Ones(gd.NxNyNz), x.cwise()*x.cwise().log() - x);
  }
  Functional grad(const Functional &, bool) const {
    assert(false);
    return 0;
  }

  void grad(const GridDescription &gd, const VectorXd &x, const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    assert(&gd); // to avoid an unused parameter error
    assert(&x); // to avoid an unused parameter error
    *outgrad += (kT*ingrad).cwise()*choose(1e-90, -207.2326583694641, x.cwise()/x + x.cwise().log() - VectorXd::Ones(gd.NxNyNz));
    if (outpgrad) *outpgrad += (kT*ingrad).cwise()*choose(1e-90, -207.2326583694641, x.cwise()/x + x.cwise().log() - VectorXd::Ones(gd.NxNyNz));
  }

  Expression printme(const Expression &) const {
    return Expression("Can't print optimized Functionals");
  }
private:
  double kT;
};

inline Functional IdealGasFast(double kT) {
  return Functional(new IdealGasFast_type(kT), "IdealGasFast");
}
