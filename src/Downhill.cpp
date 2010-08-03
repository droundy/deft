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

#include "Downhill.h"
#include <stdio.h>

bool Downhill::improve_energy(bool) {
  const VectorXd &g = grad();
  double old_energy = energy();
  *x -= nu*g;
  invalidate_cache(); // Must always remember this!!!
  if (energy() <= old_energy) {
    nu *= 1.1; // Try a bigger step next time!
  } else {
    printf("Energy went up by %g\n", energy()-old_energy);
    *x += nu*g;
    invalidate_cache();
    nu *= 0.5; // Try a smaller step next time!
  }
  return false;
}

void Downhill::print_info(int iter) const {
  Minimizer::print_info(iter);
  printf("\tnu = %g\n", nu);
}
