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

#include "ExternalPotential.h"
#include <stdio.h>
#include <math.h>

void ExternalPotential::print_summary() const {
  printf("External potential summary\n");
}

double ExternalPotential::operator()(const VectorXd &data) const {
  return Vdvolume.dot(data);
}

void ExternalPotential::grad(const VectorXd &n,
                             VectorXd *g_ptr, VectorXd *pg_ptr) const {
  *g_ptr += Vdvolume;
  if (pg_ptr) *pg_ptr += Vdvolume;
}
