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
#include <time.h>
#include <sys/utsname.h>
#include "ReciprocalGrid.h"
#include "Functionals.h"

Lattice lat(Cartesian(0,5,5), Cartesian(5,0,5), Cartesian(5,5,0));
const double R = 3, R3 = R*R*R;
int resolution = 200;
GridDescription gd(lat, resolution, resolution, resolution);
const double dr = pow(gd.dvolume, 1.0/3);

inline complex func(Reciprocal kvec) {
  double k = kvec.norm();
  double kR = k*R;
  double kdr = k*dr;
  if (kR > 1e-3) {
    return exp(-spreading*kdr*kdr)*(4*M_PI)*(sin(kR) - kR*cos(kR))/(k*k*k);
  } else {
    const double kR2 = kR*kR;
    // The following is a simple power series expansion to the above
    // function, to handle the case as k approaches zero with greater
    // accuracy (and efficiency).  I evaluate the smaller elements
    // first in the hope of reducing roundoff error (but this is not
    // yet tested).
    return (4*M_PI/3)*R3*(kR2*kR2*kR2*(-1.0/15120) + kR2*kR2*(1.0/280) + kR2*(-1.0/10) + 1 );
  }
}

clock_t start = clock();
int retval = 0;
double fft_time = 1;
double last_time;

void print_time(const char *name, double time_limit = 0) {
  clock_t stop = clock();
  last_time = (stop - double(start))/CLOCKS_PER_SEC;
  printf("%s took %g seconds which is %.2g relative to FFT.\n", name, last_time, last_time/fft_time);
  if (time_limit > 0 && last_time/fft_time > time_limit) {
    printf("FAIL: took too long!\n");
    retval++;
  }
  start = stop;
}

int main(int, char **argv) {
  Grid foo(gd);
  ReciprocalGrid recip(gd);


  print_time("Getting ready");
  recip = foo.fft();
  print_time("FFT");
  fft_time = last_time;
  recip.cwise() *= step(gd, R);
  print_time("recip.cwise() *= step(gd,R)", 2.0);
  recip = foo.fft();
  print_time("FFT", 1.5);
  for (int i=0;i<recip.rows(); i++) {
    const int z = i % gd.NzOver2;
    const int n = (i-z)/gd.NzOver2;
    const int y = n % gd.Ny;
    const int x = (n-y)/gd.Ny;
    const RelativeReciprocal rvec((x>gd.Nx/2) ? x - gd.Nx : x,
                                  (y>gd.Ny/2) ? y - gd.Ny : y,
                                  z);
    // FIXME: it seems that brillouinZone is broken...  :(
    //return func(fineLat.brillouinZone(Lat.toReciprocal(rvec)));
    recip[i] *= func(gd.Lat.toReciprocal(rvec));
  }
  print_time("recip *= step(gd,R) done by loop", 2.0);
  recip = foo.fft();
  print_time("FFT", 1.5);
  recip.MultiplyBy(func);
  print_time("recip.MultiplyBy(func)", 2.0);
  foo = StepConvolve(R)(foo);
  print_time("StepConvolve(R)(foo)", 5.0);

  if (retval == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], retval);
  return retval;
}
