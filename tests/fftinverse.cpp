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
#include "Grid.h"
#include "ReciprocalGrid.h"

double gaussian(Cartesian r) {
  const Cartesian center(0, 0, 0);
  const Cartesian off(0.5, 0, 0);
  const Cartesian off2(-0.5, 0, 0);
  const Cartesian dr(r - center), dr2(r-off), dr3(r-off2);
  return exp(-16*(dr*dr)) - 1.5*exp(-50*(dr*dr))
    - 0.5*exp(-60*(dr2*dr2)) - 0.5*exp(-60*(dr3*dr3));
}

int main(int, char *argv[]) {
  printf("Working on %s\n", argv[0]);
  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  int resolution = 10;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid foo(gd);
  foo.Set(gaussian);
  foo += 0.1*(-30*r2(gd)).cwise().exp();
  Grid foo2(foo.fft().ifft());
  int errorcode = 0;
  for (int x=0; x<resolution; x++)
    for (int y=0; y<resolution; y++)
      for (int z=0; z<resolution; z++) {
        double e = foo(x,y,z) - foo2(x,y,z);
        if (fabs(e) > 1e-15) {
          printf("Error of %g\n", e);
          errorcode += 1;
        }
      }
  return errorcode;
}
