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

#include "Functional.h"
#include "handymath.h"
#include <math.h>

// A simple integrating functional...

class Integrate : public FunctionalInterface {
public:
  Integrate(const FieldFunctional &ff) : f(ff) {};

  double energy(double data) const {
    return f(data);
  }
  double grad(double data) const {
    return f.grad(data);
  }
  double energy(const GridDescription &gd, const VectorXd &data) const {
    // This does some extra work to save the energies of each term in
    // the sum.
    VectorXd fdata(f.justMe(gd, data));
    double e = gd.dvolume*fdata.sum();
    f.last_energy = e;
    FieldFunctional *nxt = f.next();
    while (nxt) {
      fdata += nxt->justMe(gd, data);
      double etot = gd.dvolume*fdata.sum();
      nxt->last_energy = etot - e;
      e = etot;
      nxt = nxt->next();
    }
    return e;
  }
  void grad(const GridDescription &gd, const VectorXd &x, VectorXd *g, VectorXd *pgrad) const {
    f.grad(gd, x, gd.dvolume*VectorXd::Ones(gd.NxNyNz), g, pgrad);
  }
private:
  const FieldFunctional f;
};

Functional integrate(const FieldFunctional &f) {
  return Functional(new Integrate(f));
}
