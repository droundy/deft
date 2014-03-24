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

#include <stdio.h>
#include <sys/stat.h>
#include "Grid.h"
#include "Functionals.h"

Functional gSigmaA_automagic(double R);
Functional gSigmaA_by_hand(double R);

double my_density(Cartesian r) {
  const Cartesian center(1, 0, 0);
  const Cartesian dr(r - center);
  const double width = 1;
  //if (drmag > 2*width) return 0;
  const double max_density = 1/(4*M_PI/3.0);
  const double min_density = 1e-8*max_density;
  return max_density*exp(-(dr*dr)/(2*width*width)) + min_density;
}

int main(int, char **argv) {
  printf("hello world\n");
  Lattice lat(Cartesian(0.1,0,0), Cartesian(0,0.1,0), Cartesian(0,0,10));
  printf("creating grid description\n");
  GridDescription gd(lat, 0.01);
  Grid n(gd);
  n.Set(my_density);
  Grid gsigauto(gd, gSigmaA_automagic(1)(1, n));
  Grid gsighand(gd, gSigmaA_by_hand(1)(1, n));
  for (int i=0;i<gd.Nz;i++) {
    Cartesian here = gd.fineLat.toCartesian(Relative(0,0,i));
    printf("%g\t%g\t%g\t%g\t%g\n", here[2], n(here),
           gsigauto(here), gsighand(here), gsigauto(here) - gsighand(here));
  }
}
