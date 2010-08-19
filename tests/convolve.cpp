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
#include "Functionals.h"

const double themax = 4000;

double gaussian(Cartesian r) {
  const Cartesian center(0, 0, 0);
  const Cartesian dr(r - center);
  return themax*exp(-50000*(dr*dr));
}

int main(int, char **argv) {
  Lattice lat(Cartesian(0,5,5), Cartesian(5,0,5), Cartesian(5,5,0));
  Cartesian plotcorner(-5, -5, 0), plotx(10,0,0), ploty(0,10,0);
  int resolution = 100;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid foo(gd), bar(gd);
  printf("Running Set(gaussian)...\n");
  foo.Set(gaussian);
  printf("Original integrates to %.15g\n", integrate(foo));
  printf("Original Maximum is %g\n", foo.maxCoeff());
  int retval = 0;
  if (fabs(integrate(foo)-1) > 1e-6) {
    printf("Error in original is too large:  %g\n", integrate(foo)-1);
    retval++;
  }
  if (fabs(foo.maxCoeff()-themax) > 1e-6) {
    printf("Max of original is wrong:  %.15g\n", foo.maxCoeff()-themax);
    retval++;
  }
  foo.epsNativeSlice("unblurred.eps", plotx, ploty, plotcorner);

  printf("Running Gaussian(10)...\n");
  bar = Gaussian(10)(foo);
  printf("Gaussian(10) integrates to %.15g\n", integrate(bar));
  if (fabs(integrate(bar)-1) > 1e-6) {
    printf("Error in Gaussian(10) is too large:  %g\n", integrate(bar)-1);
    retval++;
  }
  bar.epsNativeSlice("gaussian-width-10.eps", plotx, ploty, plotcorner);

  printf("Gaussian(1) integrates to %g\n", integrate(bar));
  printf("Running Gaussian(1)...\n");
  bar = Gaussian(1)(foo);
  if (fabs(integrate(bar)-1) > 1e-6) {
    printf("Error in Gaussian(1) is too large:  %g\n", integrate(bar)-1);
    retval++;
  }
  bar.epsNativeSlice("gaussian-width-1.eps", plotx, ploty, plotcorner);

  printf("Running StepConvolve(1)...\n");
  bar = StepConvolve(1)(foo);
  printf("StepConvolve(1) integrates to %.15g\n", integrate(bar));
  printf("StepConvolve(1) Maximum is %g\n", bar.maxCoeff());
  if (fabs((bar.maxCoeff()-themax)/themax) > 1e-6) {
    printf("Max of StepConvolve(1) is wrong:  %g\n", (bar.maxCoeff()-themax)/themax);
    //retval++;
  }
  const double fourpiover3 = 4*M_PI/3;
  if (fabs((integrate(bar)-fourpiover3)/fourpiover3) > 1e-6) {
    printf("Integral of StepConvolve(1) is wrong:  %g\n", (integrate(bar)-fourpiover3)/fourpiover3);
    retval++;
  }
  bar.epsNativeSlice("step-1.eps", plotx, ploty, plotcorner);

  printf("Running StepConvolve(2)...\n");
  bar = StepConvolve(2)(foo);
  printf("StepConvolve(2) Maximum is %g\n", bar.maxCoeff());
  printf("StepConvolve(2) integrates to %.15g\n", integrate(bar));
  if (fabs((integrate(bar)-fourpiover3*8)/fourpiover3/8) > 1e-6) {
    printf("Integral of StepConvolve(2) is wrong:  %g\n", (integrate(bar)-fourpiover3*8)/fourpiover3/8);
    retval++;
  }
  if (fabs((bar.maxCoeff()-themax)/themax) > 1e-6) {
    printf("Max of StepConvolve(2) is wrong:  %g\n", (bar.maxCoeff()-themax)/themax);
    //retval++;
  }
  bar.epsNativeSlice("step-2.eps", plotx, ploty, plotcorner);

  printf("Running StepConvolve(3)...\n");
  bar = StepConvolve(3)(foo);
  printf("StepConvolve(3) Maximum is %g\n", bar.maxCoeff());
  printf("StepConvolve(3) integrates to %.15g\n", integrate(bar));
  if (fabs((integrate(bar)-fourpiover3*27)/fourpiover3/27) > 1e-6) {
    printf("Integral of StepConvolve(2) is wrong:  %g\n", (integrate(bar)-fourpiover3*27)/fourpiover3/27);
    retval++;
  }
  if (fabs((bar.maxCoeff()-themax)/themax) > 1e-6) {
    printf("Max of StepConvolve(3) is wrong:  %g\n", (bar.maxCoeff()-themax)/themax);
    //retval++;
  }
  bar.epsNativeSlice("step-3.eps", plotx, ploty, plotcorner);

  if (retval == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], retval);
  return retval;
}
