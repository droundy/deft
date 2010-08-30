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

double gaussian(Cartesian r) {
  const Cartesian center(1, 0, 0);
  const Cartesian dr(r - center);
  // This is smooth so that we won't suffer from oscillations due to
  // aliasing effects.  We might also evaluate the convolution
  // function in real space and do a DFT on it to get the convolution
  // in fourier space.  This would avoid the oscillations, but would
  // require that we store the convolution function on a grid.
  return 4000*exp(-50*(dr*dr));
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
  retval += integrate(Gaussian(1)).run_finite_difference_test("integrate Gaussian(1)", foo);
  retval += integrate(StepConvolve(1)).run_finite_difference_test("integrate StepConvolve(1)", foo);
  retval += integrate(StepConvolve(3)).run_finite_difference_test("integrate StepConvolve(3)", foo);
  retval += integrate(ShellConvolve(1)).run_finite_difference_test("integrate ShellConvolve(1)", foo);
  retval += integrate(ShellConvolve(3)).run_finite_difference_test("integrate ShellConvolve(3)", foo);
  FieldFunctional ysh = yShellConvolve(1), sh = ShellConvolve(1);
  retval += integrate(sqr(ysh)).run_finite_difference_test("integrate ysh^2", foo);
  retval += integrate(ysh*ysh).run_finite_difference_test("integrate ysh*ysh", foo);
  retval += integrate(ysh*ysh + sh).run_finite_difference_test("integrate ysh*ysh + sh", foo);
  retval += integrate(-1*ysh*ysh).run_finite_difference_test("integrate -1*ysh*ysh", foo);
  retval += (-1*integrate(ysh*ysh)).run_finite_difference_test("-1* integrate ysh*ysh", foo);

  FieldFunctional zxsh = zxShellConvolve(3);
  retval += integrate(sqr(sqr(zxsh))).run_finite_difference_test("integrate zxsh^2", foo);
  {
    Grid foo2(foo);
    foo2.cwise() -= 1; // so we don't divide by zero...
    retval += integrate(sqr(zxsh)/sh).run_finite_difference_test("integrate zxsh^2/sh", foo2);
  }

  if (retval == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], retval);
  return retval;
}
