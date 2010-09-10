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

class ExternalPotentialType : public FieldFunctionalInterface {
public:
  explicit ExternalPotentialType(const VectorXd &Vin) : V(Vin) {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    return V.cwise()*data;
  }
  double transform(double) const {
    return 0;
  }
  double grad(double) const {
    return 0;
  }

  void grad(const GridDescription &, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += V.cwise()*ingrad;
    if (outpgrad) *outpgrad += V.cwise()*ingrad;
  }
private:
  VectorXd V; // the potential times the differential volume element dV.
};

FieldFunctional ExternalPotential(const VectorXd &V) {
  return FieldFunctional(new ExternalPotentialType(V), "external potential");
}
