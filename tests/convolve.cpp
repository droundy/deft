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

Functional TensorDensityXX(double R);
Functional VectorDensityX(double R);
Functional n2Density(double R);

double gaussian(Cartesian r) {
  const Cartesian center(1, 0, 0);
  const Cartesian dr(r - center);
  // This is smooth so that we won't suffer from oscillations due to
  // aliasing effects.  We might also evaluate the convolution
  // function in real space and do a DFT on it to get the convolution
  // in fourier space.  This would avoid the oscillations, but would
  // require that we store the convolution function on a grid.
  //return 4000*exp(-50*(dr*dr));

  // Actually, it's sharp, because we want to make sure the
  // oscillations from aliasing are under controll.
  return 4000*exp(-50*1000.0*(dr*dr));
}

int main(int, char **argv) {
  const double kT = 1e-3;
  Lattice lat(Cartesian(0,5,5), Cartesian(5,0,5), Cartesian(5,5,0));
  Cartesian plotcorner(-5, -5, 0), plotx(10,0,0), ploty(0,10,0);
  int resolution = 200;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid foo(gd), bar(gd), ref(gd);
  printf("Running Set(gaussian)...\n");
  foo.Set(gaussian);
  const double integrate_foo = Identity().integral(kT, foo);
  printf("Original integrates to %.15g\n", integrate_foo);
  printf("Original Maximum is %g\n", foo.maxCoeff());
  int retval = 0;
  if (integrate_foo < 0) {
    printf("Integral of original is negative (which may throw off tests):  %g\n", integrate_foo);
    retval++;
  }
  //mkdir("tests/vis", 0777);
  //foo.epsNativeSlice("tests/vis/unblurred.eps", plotx, ploty, plotcorner);

  printf("Running Gaussian(10)...\n");
  bar = Gaussian(10)(kT, foo);
  printf("Gaussian(10) integrates to %.15g\n", Identity().integral(kT, bar));
  if (fabs(Identity().integral(kT, bar)/integrate_foo-1) > 1e-6) {
    printf("Error in Gaussian(10) is too large:  %g\n", Identity().integral(kT, bar)/integrate_foo-1);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/gaussian-width-10.eps", plotx, ploty, plotcorner);

  printf("Gaussian(1) integrates to %g\n", Identity().integral(kT, bar));
  printf("Running Gaussian(1)...\n");
  bar = Gaussian(1)(kT, foo);
  if (fabs(Identity().integral(kT, bar)/integrate_foo-1) > 1e-6) {
    printf("Error in Gaussian(1) is too large:  %g\n", Identity().integral(kT, bar)/integrate_foo-1);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/gaussian-width-1.eps", plotx, ploty, plotcorner);

  printf("Running StepConvolve(1)...\n");
  bar = StepConvolve(1)(kT, foo);
  printf("StepConvolve(1) integrates to %.15g\n", Identity().integral(kT, bar));
  printf("StepConvolve(1) Maximum is %g\n", bar.maxCoeff());
  if (fabs(bar.maxCoeff()/integrate_foo-1) > 1e-6) {
    printf("FAIL: Max of StepConvolve(1) is wrong:  %g\n", bar.maxCoeff()/integrate_foo - 1);
    retval++;
  }
  printf("StepConvolve(1) Minimum is %g\n", bar.minCoeff());
  if (bar.minCoeff()/integrate_foo < -1e-9) {
    printf("Min of StepConvolve(1) is wrong:  %g\n", bar.minCoeff()/integrate_foo);
    retval++;
  }
  const double fourpiover3 = 4*M_PI/3;
  if (fabs((Identity().integral(kT, bar)/integrate_foo-fourpiover3)/fourpiover3) > 1e-6) {
    printf("Integral of StepConvolve(1) is wrong:  %g\n",
           (Identity().integral(kT, bar)/integrate_foo-fourpiover3)/fourpiover3);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/step-1.eps", plotx, ploty, plotcorner);

  printf("Running StepConvolve(2)...\n");
  bar = StepConvolve(2)(kT, foo);
  printf("StepConvolve(2) Maximum is %g\n", bar.maxCoeff());
  printf("StepConvolve(2) integrates to %.15g\n", Identity().integral(kT, bar));
  if (fabs((Identity().integral(kT, bar)/integrate_foo-fourpiover3*8)/fourpiover3/8) > 1e-6) {
    printf("Integral of StepConvolve(2) is wrong:  %g\n", (Identity().integral(kT, bar)/integrate_foo-fourpiover3*8)/fourpiover3/8);
    retval++;
  }
  if (fabs((bar.maxCoeff()-integrate_foo)/integrate_foo) > 1e-6) {
    printf("Max of StepConvolve(2) is wrong:  %g\n", (bar.maxCoeff()-integrate_foo)/integrate_foo);
    retval++;
  }
  printf("StepConvolve(2) Minimum is %g\n", bar.minCoeff());
  if (bar.minCoeff()/integrate_foo < -1e-9) {
    printf("Min of StepConvolve(2) is wrong:  %g\n", bar.minCoeff()/integrate_foo);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/step-2.eps", plotx, ploty, plotcorner);

  printf("Running StepConvolve(3)...\n");
  bar = StepConvolve(3)(kT, foo);
  printf("StepConvolve(3) Maximum is %g\n", bar.maxCoeff());
  printf("StepConvolve(3) integrates to %.15g\n", Identity().integral(kT, bar));
  if (fabs((Identity().integral(kT, bar)/integrate_foo-fourpiover3*27)/fourpiover3/27) > 1e-6) {
    printf("Integral of StepConvolve(3) is wrong:  %g\n",
           (Identity().integral(kT, bar)/integrate_foo-fourpiover3*27)/fourpiover3/27);
    retval++;
  }
  if (fabs((bar.maxCoeff()-integrate_foo)/integrate_foo) > 1e-6) {
    printf("Max of StepConvolve(3) is wrong:  %g\n", (bar.maxCoeff()-integrate_foo)/integrate_foo);
    retval++;
  }
  printf("StepConvolve(3) Minimum is %g\n", bar.minCoeff());
  if (fabs(bar.minCoeff()/integrate_foo) < -1e-9) {
    printf("Min of StepConvolve(3) is wrong:  %g\n", bar.minCoeff()/integrate_foo);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/step-3.eps", plotx, ploty, plotcorner);


  const double fourpi = 4*M_PI;
  printf("Running ShellConvolve(1)...\n");
  bar = ShellConvolve(1)(kT, foo);
  printf("ShellConvolve(1) integrates to %.15g\n", Identity().integral(kT, bar));
  printf("ShellConvolve(1) Maximum is %g\n", bar.maxCoeff());
  if (fabs((Identity().integral(kT, bar)/integrate_foo-fourpi)/fourpi) > 1e-6) {
    printf("Integral of ShellConvolve(1) is wrong:  %g\n",
           (Identity().integral(kT, bar)/integrate_foo-fourpi)/fourpi);
    retval++;
  }
  printf("ShellConvolve(1) Minimum is %g\n", bar.minCoeff());
  if (bar.minCoeff()/integrate_foo < -1e-9) {
    printf("Min of ShellConvolve(1) is wrong:  %g\n", bar.minCoeff()/integrate_foo);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/shell-1.eps", plotx, ploty, plotcorner);


  printf("Running ShellPrimeConvolve(1)...\n");
  bar = ShellPrimeConvolve(1)(kT, foo);
  printf("ShellPrimeConvolve(1) integrates to %.15g\n", Identity().integral(kT, bar));
  printf("ShellPrimeConvolve(1) Maximum is %g\n", bar.maxCoeff());
  //bar.epsNativeSlice("tests/vis/shellPrime-1.eps", plotx, ploty, plotcorner);

  printf("Running ShellConvolve(3)...\n");
  bar = ShellConvolve(3)(kT, foo);
  printf("ShellConvolve(3) Maximum is %g\n", bar.maxCoeff());
  printf("ShellConvolve(3) integrates to %.15g\n", Identity().integral(kT, bar));
  if (fabs((Identity().integral(kT, bar)/integrate_foo-fourpi*9)/fourpi/9) > 1e-6) {
    printf("Integral of ShellConvolve(3) is wrong:  %g\n",
           (Identity().integral(kT, bar)/integrate_foo-fourpi*9)/fourpi/9);
    retval++;
  }
  printf("ShellConvolve(3) Minimum is %g\n", bar.minCoeff());
  if (bar.minCoeff()/integrate_foo < -1e-9) {
    printf("Min of ShellConvolve(3) is wrong:  %g\n", bar.minCoeff()/integrate_foo);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/shell-3.eps", plotx, ploty, plotcorner);

  printf("Running yShellConvolve(1)...\n");
  bar = yShellConvolve(1)(kT, foo);
  printf("yShellConvolve(1) integrates to %.15g\n", Identity().integral(kT, bar));
  printf("yShellConvolve(1) Maximum is %g\n", bar.maxCoeff());
  printf("yShellConvolve(1) Minimum is %g\n", bar.minCoeff());
  if (fabs(Identity().integral(kT, bar)/integrate_foo) > 1e-6) {
    printf("Integral of yShellConvolve(1) is wrong:  %g\n",
           Identity().integral(kT, bar)/integrate_foo-fourpi);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/y-shell-1.eps", plotx, ploty, plotcorner);
  Functional ysh = yShellConvolve(1), sh = ShellConvolve(1);
  {
    Grid scalarsh(gd, ShellConvolve(1)(kT, foo));
    if (bar.maxCoeff() > scalarsh.maxCoeff()) {
      printf("FAIL: max of vector shell greater than scalar shell!\n");
      retval++;
    }
    if (-bar.minCoeff() > scalarsh.maxCoeff()) {
      printf("FAIL: min of vector shell greater in magnitude than scalar shell!\n");
      retval++;
    }
  }

  printf("Running xShellConvolve(3)...\n");
  ref = ShellConvolve(3)(kT, foo);
  bar = xShellConvolve(3)(kT, foo);
  printf("xShellConvolve(3) Maximum is %g\n", bar.maxCoeff());
  printf("xShellConvolve(3) integrates to %.15g\n", Identity().integral(kT, bar));
  if (fabs(Identity().integral(kT, bar)/integrate_foo) > 1e-6) {
    printf("Integral of xShellConvolve(3) is wrong:  %g\n",
           Identity().integral(kT, bar)/integrate_foo);
    retval++;
  }
  double mymax = ref.maxCoeff();
  if (fabs(bar.maxCoeff() - mymax)/mymax > 2e-3) {
    printf("problem with xShellConvolve max: %g != %g (error of %g)\n",
           bar.maxCoeff(), mymax, (bar.maxCoeff() - mymax)/mymax);
    retval++;
    exit(1);
  }
  if (fabs((bar.minCoeff() + mymax)/mymax) > 2e-3) {
    printf("problem with xShellConvolve min: %g != minus %g (error of %g)\n",
           bar.minCoeff(), mymax, (bar.minCoeff() + mymax)/mymax);
    retval++;
    exit(1);
  }
  //bar.epsNativeSlice("tests/vis/x-shell-3.eps", plotx, ploty, plotcorner);
  for (int i=0;i<gd.NxNyNz;i++) {
    if ((fabs(bar[i]) - fabs(ref[i]))/mymax > 1e-11) {
      printf("FAIL: x shell is bigger in magnitude than shell by %g out of %g compared with %g!\n",
             fabs(bar[i]) - fabs(ref[i]), fabs(ref[i]), mymax);
      retval++;
    }
  }

  printf("Running VectorDensityX(3)...\n");
  ref = n2Density(3)(kT, foo);
  bar = VectorDensityX(3)(kT, foo);
  printf("VectorDensityX(3) Maximum is %g\n", bar.maxCoeff());
  printf("VectorDensityX(3) integrates to %.15g\n", Identity().integral(kT, bar));
  if (fabs(Identity().integral(kT, bar)/integrate_foo) > 1e-6) {
    printf("Integral of VectorDensityX(3) is wrong:  %g\n",
           Identity().integral(kT, bar)/integrate_foo);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/x-vector-3.eps", plotx, ploty, plotcorner);
  mymax = ref.maxCoeff();
  if (fabs(bar.maxCoeff() - mymax)/mymax > 2e-3) {
    printf("problem with VectorDensityX max: %g != %g (error of %g)\n",
           bar.maxCoeff(), mymax, (bar.maxCoeff() - mymax)/mymax);
    retval++;
    exit(1);
  }
  if (fabs((bar.minCoeff() + mymax)/mymax) > 2e-3) {
    printf("problem with VectorDensityX min: %g != minus %g (error of %g)\n",
           bar.minCoeff(), mymax, (bar.minCoeff() + mymax)/mymax);
    retval++;
    exit(1);
  }
  for (int i=0;i<gd.NxNyNz;i++) {
    if ((fabs(bar[i]) - fabs(ref[i]))/mymax > 1e-11) {
      printf("FAIL: x shell is bigger in magnitude than shell by %g out of %g compared with %g!\n",
             fabs(bar[i]) - fabs(ref[i]), fabs(ref[i]), mymax);
      retval++;
    }
  }

  printf("Running xyShellConvolve(1)...\n");
  ref = xShellConvolve(1)(kT, foo);
  bar = xyShellConvolve(1)(kT, foo);
  printf("xyShellConvolve(1) integrates to %.15g\n", Identity().integral(kT, bar));
  printf("xyShellConvolve(1) Maximum is %g\n", bar.maxCoeff());
  if (fabs(Identity().integral(kT, bar)/integrate_foo) > 1e-6) {
    printf("Integral of xyShellConvolve(1) is wrong:  %g\n",
           Identity().integral(kT, bar)/integrate_foo-fourpi);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/xy-shell-1.eps", plotx, ploty, plotcorner);
  mymax = ref.maxCoeff();
  for (int i=0;i<gd.NxNyNz;i++) {
    if ((fabs(bar[i]) - fabs(ref[i]))/mymax > 1e-11) {
      printf("FAIL: xy shell is bigger in magnitude than x shell by %g out of %g!\n",
             fabs(bar[i]) - fabs(ref[i]), fabs(ref[i]));
      retval++;
    }
  }

  printf("Running xxShellConvolve(2)...\n");
  ref = xShellConvolve(2)(kT, foo).cwise().abs() + 1./3*ShellConvolve(2)(kT, foo);
  bar = xxShellConvolve(2)(kT, foo);
  printf("xxShellConvolve(2) integrates to %.15g\n", Identity().integral(kT, bar));
  printf("xxShellConvolve(2) Maximum is %g\n", bar.maxCoeff());
  if (fabs(Identity().integral(kT, bar)/integrate_foo) > 1e-14) {
    printf("FAIL: Integral of xxShellConvolve(2) is wrong:  %g\n",
           Identity().integral(kT, bar)/integrate_foo);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/xx-shell-2.eps", plotx, ploty, plotcorner);
  //bar.epsNativeSlice("tests/vis/xx-shell-2.eps", Cartesian(0,10,0), Cartesian(0,0,10), Cartesian(0,0,0));
  mymax = ref.maxCoeff();
  for (int i=0;i<gd.NxNyNz;i++) {
    if ((fabs(bar[i]) - fabs(ref[i]))/mymax > 1e-12) {
      printf("FAIL: xx shell is bigger in magnitude than x shell by %g out of %g compared with %g!\n",
             fabs(bar[i]) - fabs(ref[i]), fabs(ref[i]), mymax);
      retval++;
    }
  }

  printf("Running TensorDensityXX(2)...\n");
  ref = VectorDensityX(2)(kT, foo).cwise().abs() + 1./3*n2Density(2)(kT, foo);
  bar = TensorDensityXX(2)(kT, foo);
  printf("TensorDensityXX(2) integrates to %.15g\n", Identity().integral(kT, bar));
  printf("TensorDensityXX(2) Maximum is %g\n", bar.maxCoeff());
  if (fabs(Identity().integral(kT, bar)/integrate_foo) > 1e-14) {
    printf("FAIL: Integral of xxShellConvolve(2) is wrong:  %g\n",
           Identity().integral(kT, bar)/integrate_foo);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/xx-tensor-2.eps", plotx, ploty, plotcorner);
  mymax = ref.maxCoeff();
  double shellmax = n2Density(3)(kT, foo).maxCoeff();
  if (fabs(bar.maxCoeff() - 2*shellmax/3)/(2*shellmax/3) > 1e-2) {
    printf("problem with TensorDensityXX max: %g != %g (error of %g)\n",
           bar.maxCoeff(), 2*shellmax/3, (bar.maxCoeff() - 2*shellmax/3)/(2*shellmax/3));
    retval++;
    exit(1);
  }
  for (int i=0;i<gd.NxNyNz;i++) {
    if ((fabs(bar[i]) - fabs(ref[i]))/mymax > 1e-12) {
      printf("FAIL: xx shell is bigger in magnitude than x shell by %g out of %g compared with %g!\n",
             fabs(bar[i]) - fabs(ref[i]), fabs(ref[i]), mymax);
      retval++;
    }
  }
  printf("Running zxShellConvolve(3)...\n");
  Functional zxsh = zxShellConvolve(3);
  bar = zxsh(kT, foo);
  printf("zxShellConvolve(3) Maximum is %g\n", bar.maxCoeff());
  printf("zxShellConvolve(3) integrates to %.15g\n", Identity().integral(kT, bar));
  if (fabs(Identity().integral(kT, bar)/integrate_foo) > 1e-6) {
    printf("Integral of zxShellConvolve(3) is wrong:  %g\n",
           Identity().integral(kT, bar)/integrate_foo);
    retval++;
  }
  //bar.epsNativeSlice("tests/vis/zx-shell-3.eps", plotx, ploty, Cartesian(-5,-5,0.5));

  if (retval == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], retval);
  return retval;
}
