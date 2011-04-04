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
#include "CallMe.h"

// The following is a utility function used to create debug prints on
// functionals.
class CallMeClass : public FunctionalInterface {
public:
  CallMeClass(const Functional &myf, const std::string mypattern, const std::string myargs)
    : f(myf), pattern(mypattern), args(myargs) {};
  bool append_to_name(const std::string) {
    return false;
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &kT, const VectorXd &data) const {
    return f(gd, kT, data);
  }
  double transform(double kT, double n) const {
    return f(kT, n);
  }
  double derive(double kT, double n) const {
    return f.derive(kT, n);
  }
  Expression derive_homogeneous(const Expression &kT, const Expression &x) const {
    return f.derive_homogeneous(kT, x);
  }
  double d_by_dT(double kT, double n) const {
    return f.d_by_dT(kT, n);
  }
  Functional grad(const Functional &ingrad, const Functional &n, bool) const {
    if (ingrad.I_am_homogeneous() || f.I_am_local())
      // the following makes the assumption that ingrad is uniform!
      return ingrad*CallMe(f.grad(Functional(1), Identity(), false), pattern + "Grad", args)(n);
    else
      return f.grad(ingrad, n, false).set_name_stdstring(pattern + "Grad_slow");
  }
  Functional grad_T(const Functional &ingradT) const {
    return ingradT*CallMe(f.grad_T(Functional(1)), pattern + "_by_dT", args);
  }
  void grad(const GridDescription &gd, const VectorXd &kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    f.grad(gd, kT, data, ingrad, outgrad, outpgrad);
  }
  Expression printme(const Expression &kT, const Expression &x) const {
    if (x.typeIs("double")) return funexpr((pattern + args).c_str(), kT, x).set_type("double");
    return funexpr((pattern + args).c_str(), Expression("gd"), kT, x);
  }
private:
  Functional f;
  const std::string pattern, args;
};



Functional CallMe(Functional f, const std::string namepattern, const std::string args) {
  return Functional(new CallMeClass(f, namepattern, args)).set_name_stdstring(namepattern+"_value");
}
