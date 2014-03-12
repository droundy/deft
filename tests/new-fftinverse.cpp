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
#include <math.h>
#include "new/Vector.h"

int main(int, char *argv[]) {
  printf("Working on %s\n", argv[0]);
  const int Nx = 100, Ny = 98, Nz = 96;
  const int NxNyNz = Nx*Ny*Nz;
  Vector rs(NxNyNz);

  const double volume = 397;
  const double dV = volume / NxNyNz;
  const double dx = pow(dV, 1.0/3);
  for (int i=0; i<NxNyNz; i++) {
    rs[i] = exp(-(i*dx)*(i*dx)/1000);
  }
  ComplexVector ks = fft(Nx, Ny, Nz, dV, rs);
  Vector rs_again = ifft(Nx, Ny, Nz, dV, ks);
  int errorcode = 0;
  for (int i=0;i<NxNyNz;i++) {
    double e = rs[i] - rs_again[i];
    if (fabs(e) > 1e-15) {
      printf("Error of %g\n", e);
      errorcode += 1;
    }
  }
  return errorcode;
}
