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
#include "RenameAppending.h"

// The following is a utility function used to create debug prints on
// functionals.
class RenameAppendingClass : public FunctionalInterface {
public:
  RenameAppendingClass(const Functional &myf, const char *myname)
    : f(myf), name(myname) {
    printf("I am in RenameAppendingClass for %s\n", myname);
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
  Functional grad(const Functional &ingrad, const Functional &n, bool ispgrad) const {
    printf("I'm in RenameAppendingClass grad\n");
    return f.grad(ingrad, n, ispgrad);
  }
  Functional grad_T(const Functional &ingradT) const {
    printf("I'm in RenameAppendingClass grad_T\n");
    return f.grad_T(ingradT);
  }
  void grad(const GridDescription &gd, const VectorXd &kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    f.grad(gd, kT, data, ingrad, outgrad, outpgrad);
  }
  Expression printme(const Expression &kT, const Expression &x) const {
    return f.printme(kT, x).append_to_alias(name);
  }
private:
  Functional f;
  const char *name;
};


Functional RenameAppending(Functional f, const char *name) {
  return Functional(new RenameAppendingClass(f, name));
}
