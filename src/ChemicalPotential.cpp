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

class ChemicalPotentialType : public FunctionalInterface {
public:
  ChemicalPotentialType(const GridDescription &g, double chemical_potential)
    : gd(g), mu(chemical_potential) {}

  double energy(const VectorXd &data) const;
  double energy(double data) const;

  void grad(const VectorXd &data,
            VectorXd *, VectorXd *pgrad = 0) const;
private:
  GridDescription gd;
  double mu; // the chemical potential
};

double ChemicalPotentialType::energy(const VectorXd &data) const {
  double Ntot = 0;
  for (int i=0; i < gd.NxNyNz; i++) Ntot += data[i];
  Ntot *= gd.Lat.volume()/gd.NxNyNz;
  return mu*Ntot;
}

double ChemicalPotentialType::energy(double n) const {
  return mu*n;
}

void ChemicalPotentialType::grad(const VectorXd &, VectorXd *g_ptr, VectorXd *pg_ptr) const {
  VectorXd &g = *g_ptr;

  const double mudV = mu*gd.Lat.volume()/gd.NxNyNz;
  for (int i=0; i < gd.NxNyNz; i++) g[i] += mudV;
  if (pg_ptr)
    for (int i=0; i < gd.NxNyNz; i++) (*pg_ptr)[i] += mudV;
}

Functional ChemicalPotential(const GridDescription &g, double chemical_potential) {
  return Functional(new ChemicalPotentialType(g, chemical_potential));
}
