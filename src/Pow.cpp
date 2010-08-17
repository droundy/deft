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

class PowType : public FieldFunctionalInterface {
public:
  PowType(int nn) : n(nn) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    VectorXd out(data);
    for (int i=0; i<gd.NxNyNz; i++)
      for (int p=1; p < n; p++)
        out[i] *= data[i];
    return out;
  }
  double transform(double x) const {
    double v = x;
    for (int p=1; p < n; p++) v *= x;
    return v;
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
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
        (*outgrad)[i] += foo;
      }
  }
private:
  int n;
};

FieldFunctional Pow(int n) {
  assert(n > 0);
  return FieldFunctional(new PowType(n));
}
