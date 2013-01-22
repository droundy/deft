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
#include "ReciprocalGrid.h"

double gaussian(Cartesian r) {
  const Cartesian center(0.25, 0.25, .25);
  const Cartesian dr(r - center);
  return -10*exp(-500*(dr*dr));
}

int main() {
  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  //Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,1));
  Cartesian middle(0.5,0.5,0.5);
  Cartesian plotcorner(-0.5, -0.5, 0.5);
  Relative middlerel = lat.toRelative(middle);
  Reciprocal recip(0.2,0,0);
  int resolution = 20;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid foo(gd), bar(gd);
  ReciprocalGrid rfoo(gd);
  foo.Set(gaussian);
  foo += 10*(-500*r2(gd)).cwise().exp();
  //foo += 0.5*foo.x();
  mkdir("tests/vis", 0777);
  foo.epsSlice("tests/vis/demo.eps", Cartesian(1,0,0), Cartesian(0,1,0), plotcorner, 150);
  foo.epsNativeSlice("tests/vis/native.eps", Cartesian(1,0,0), Cartesian(0,1,0),
                     plotcorner);
  foo.fft().ifft().epsNativeSlice("tests/vis/native-ffted.eps", Cartesian(1,0,0),
                                  Cartesian(0,1,0), plotcorner);

  rfoo = foo.fft();
  rfoo.cwise() *= (-0.1*g2(gd)).cwise().exp();
  foo = rfoo.ifft();
  foo.epsNativeSlice("tests/vis/native-blurred.eps", Cartesian(1,0,0), Cartesian(0,1,0),
                     plotcorner);
  rfoo = (-0.4*g2(gd)).cwise().exp();
  rfoo.ifft().epsNativeSlice("tests/vis/gaussian.eps", Cartesian(1,0,0), Cartesian(0,1,0),
                             plotcorner);

  //std::cout << "and here is the foo" << foo << std::endl;
  //std::cout << "and here is the foo" << (foo + 2*bar + foo.cwise()*bar);
  std::cout << "middle and middlerel are:\n"
            << middle << std::endl << middlerel << std::endl
            << lat.toCartesian(middlerel) << std::endl;
  std::cout << "middle dot recip: " << recip * middle << std::endl;
  return 0;
}
