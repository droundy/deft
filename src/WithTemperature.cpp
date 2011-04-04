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

class WithTemperatureClass : public FunctionalInterface {
public:
  WithTemperatureClass(const Functional &newkT, const Functional &myf)
    //: mykT(newkT), f(myf), f_grad_T_1(1) {
    : mykT(newkT), f(myf), f_grad_T_1(myf.grad_T(Functional(1))) {
    printf("I am in WithTemperatureClass constructor\n");
  }
  bool append_to_name(const std::string x) {
    mykT.append_to_name(x);
    f.append_to_name(x);
    return true;
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &kT, const VectorXd &data) const {
    return f(gd, mykT(gd, kT, data), data);
  }
  double transform(double kT, double n) const {
    return f(mykT(kT, n), n);
  }
  double derive(double kT, double n) const {
    double kTnew = mykT(kT, n);
    // d/dn = d/dn + d/dT * dTnew/dn
    return f.derive(kTnew, n) + f.d_by_dT(kTnew, n)*mykT.derive(kT, n);
  }
  Expression derive_homogeneous(const Expression &kT, const Expression &x) const {
    Expression kTnew = mykT.printme(kT, x);
    return f.derive_homogeneous(kTnew, x)
      + f_grad_T_1.printme(kTnew, x)*mykT.derive_homogeneous(kT, x);
  }
  double d_by_dT(double kT, double n) const {
    double kTnew = mykT(kT, n);
    // d/dn = d/dn + d/dT * dTnew/dn
    return f.d_by_dT(kTnew, n)*mykT.d_by_dT(kT, n);
  }
  Functional grad(const Functional &ingrad, const Functional &n, bool ispgrad) const {
    printf("calling grad in WithTemperatureClass...\n");
    return WithTemperature(mykT, f.grad(ingrad, n, ispgrad))
      + mykT.grad(f.grad_T(ingrad), n, ispgrad);
  }
  Functional grad_T(const Functional &ingradT) const {
    printf("calling grad_T in WithTemperatureClass...\n");
    return mykT.grad_T(WithTemperature(mykT, f.grad_T(ingradT)));
  }
  void grad(const GridDescription &gd, const VectorXd &kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    VectorXd kTnew = mykT(gd, kT, data);
    f.grad(gd, kTnew, data, ingrad, outgrad, outpgrad);
    // Now include gradient of any density-dependence in the effective
    // temperature:
    mykT.grad(gd, kT, data, f_grad_T_1(gd, kTnew, data).cwise()*ingrad, outgrad, outpgrad);
  }
  Expression printme(const Expression &kT, const Expression &x) const {
    return f.printme(mykT.printme(kT, x), x);
  }
private:
  Functional mykT, f, f_grad_T_1;
};

Functional WithTemperature(const Functional &newkT, const Functional &f) {
  return Functional(new WithTemperatureClass(newkT, f));
}
