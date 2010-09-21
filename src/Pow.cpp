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

class PowType : public FunctionalInterface {
public:
  PowType(int nn) : n(nn) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    switch (n) {
    case 0: return VectorXd::Ones(data.rows());
    case 1: return data;
    }
    VectorXd out(data);
    for (int i=0; i<gd.NxNyNz; i++)
      for (int p=1; p < n; p++)
        out[i] *= data[i];
    return out;
  }
  double transform(double x) const {
    switch (n) {
    case 0: return 1;
    case 1: return x;
    }
    double v = x;
    for (int p=1; p < n; p++) v *= x;
    return v;
  }
  double derive(double x) const {
    if (n < 1) return 0;
    double v = n;
    for (int p=1; p < n; p++) v *= x;
    return v;
  }
  Functional grad(const Functional &ingrad, bool) const {
    switch (n) {
    case 0: return 0;
    case 1: return ingrad;
    case 2: return 2*Identity()*ingrad;
    }
    return n*Pow(n-1)*ingrad;
  }
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    switch (n) {
    case 0: return; // zero gradient!
    case 1:
      *outgrad += ingrad;
      if (outpgrad) *outpgrad += ingrad;
      return;
    }
    for (int i=0; i<gd.NxNyNz; i++) {
      double foo = n*ingrad[i];
      for (int p=1; p < n; p++) foo *= data[i];
      (*outgrad)[i] += foo;
    }
    // FIXME: we will want to propogate preexisting preconditioning eventually...
    if (outpgrad)
      for (int i=0; i<gd.NxNyNz; i++) {
        double foo = n*ingrad[i];
        for (int p=1; p < n; p++) foo *= data[i];
        (*outpgrad)[i] += foo;
      }
  }
  Expression printme(const Expression &x) const {
    switch (n) {
    case 0: return 0;
    case 1: return x;
    case 2: return x.method("cwise").method("square");
    }
    // This is more than a little hokey...
    if (n & 1) return x*Pow(n-1).printme(x);
    return Pow(n/2).printme(x.method("cwise").method("square"));
  }
private:
  int n;
};

Functional Pow(int n) {
  assert(n > 0);
  return Functional(new PowType(n));
}
