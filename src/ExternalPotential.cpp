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

class ExternalPotentialType : public FunctionalInterface {
public:
  ExternalPotentialType(const Grid &V)
    : Vdvolume(V*V.description().Lat.volume()/V.description().NxNyNz) {}

  double energy(const VectorXd &data) const {
    return Vdvolume.dot(data);
  }
  double energy(double) const {
    return Vdvolume.sum(); // FIXME: this is broken!
  }
  void grad(const VectorXd &,VectorXd *g_ptr, VectorXd *pg_ptr = 0) const {
    *g_ptr += Vdvolume;
    if (pg_ptr) *pg_ptr += Vdvolume;
  }

  void print_summary(const char *prefix, const VectorXd &data) const {
    printf("%sExternal potential energy = %g\n", prefix, (*this)(data));
  }
private:
  VectorXd Vdvolume; // the potential times the differential volume element dV.
};

Functional ExternalPotential(const Grid &V) {
  return Functional(new ExternalPotentialType(V));
}
