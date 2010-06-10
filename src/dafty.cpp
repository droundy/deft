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
#include "Grid.h"
#include "ReciprocalGrid.h"

double gaussian(Cartesian r) {
  const Cartesian center(0.125, 0.125, 0);
  const Cartesian dr(r - center);
  return -0*exp(-500*(dr*dr));
}

int main() {
  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  Cartesian middle(0.5,0.5,0.5);
  Relative middlerel = lat.toRelative(middle);
  Reciprocal recip(0.2,0,0);
  int resolution = 20;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid foo(gd), bar(gd);
  ReciprocalGrid rfoo(gd);
  foo.Set(gaussian);
  foo += 10*(-500*foo.r2()).cwise().exp();
  //foo += 0.5*foo.x();
  foo.epsSlice("demo.eps", Cartesian(1,0,0), Cartesian(0,1,0),
               Cartesian(-.5,-.5,.5), 150);
  foo.epsNativeSlice("native.eps", Cartesian(1,0,0), Cartesian(0,1,0),
                     Cartesian(-.5,-.5,.5));
  foo.fft().ifft().epsNativeSlice("native-ffted.eps",
                                  Cartesian(1,0,0), Cartesian(0,1,0),
                                  Cartesian(-.5,-.5,.5));

  rfoo = foo.fft();
  rfoo.cwise() *= (-10.1*rfoo.g2()).cwise().exp();
  foo = rfoo.ifft();
  foo.epsNativeSlice("native-blurred.eps", Cartesian(1,0,0), Cartesian(0,1,0),
                     Cartesian(-.5,-.5,.5));
  
  //std::cout << "and here is the foo" << foo << std::endl;
  //std::cout << "and here is the foo" << (foo + 2*bar + foo.cwise()*bar);
  std::cout << "middle and middlerel are:\n"
            << middle << std::endl << middlerel << std::endl
            << lat.toCartesian(middlerel) << std::endl;
  std::cout << "middle dot recip: " << recip * middle << std::endl;
  return 0;
}
