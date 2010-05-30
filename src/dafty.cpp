// Dafty is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Dafty Authors
//
// Dafty is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with dafty; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include <stdio.h>
#include "grid.h"

int main() {
  Lattice lat(Cartesian(0,1,1), Cartesian(1,0,1), Cartesian(1,1,0));
  Cartesian middle(0.5,0.5,0.5);
  Relative middlerel = lat.toRelative(middle);
  Reciprocal recip(0.2,0,0);
  Grid foo(lat, 4, 4, 4), bar(lat, 4, 4, 4);
  std::cout << "and here is the foo" << (foo + bar);
  std::cout << "middle and middlerel are:\n"
            << middle << std::endl << middlerel << std::endl
            << lat.toCartesian(middlerel) << std::endl;
  std::cout << "middle dot recip: " << recip * middle << std::endl;
  return 0;
}
