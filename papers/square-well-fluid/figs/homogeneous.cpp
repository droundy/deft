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
#include <cassert>
#include <sys/stat.h>
#include <sys/types.h>

#include "new/HomogeneousSW_liquidFast.h"
#include "version-identifier.h"
#include "handymath.h"

// The following tells fac how to run the executable to generate a
// homogeneous.dat file.
/*--

for kT in list(range(1,11))+[100]:
    for ww in [1.3]:
        self.rule('%s %g %g' % (exe,ww,kT), [exe],
                  [os.path.dirname(exe)[:-5]+'/data/homogeneous/ww%g-kT%g.dat' % (ww,kT)])

--*/
const double epsilon = 1.0;
const double radius = 1.0;
const double sigma = 2*radius;

int main(int argc, char **argv) {
  double temp, lambda;
  if (argc != 3) {
    printf("usage: %s lambda kT\n", argv[0]);
    return 1;
  }
  printf("git version: %s\n", version_identifier());

  sscanf(argv[1], "%lg", &lambda);
  sscanf(argv[2], "%lg", &temp);

  printf("lam %g temp %g\n", lambda, temp);

  char *fname = new char[4096];
  mkdir("papers/square-well-fluid/data/homogeneous", 0777);
  sprintf(fname, "papers/square-well-fluid/data/homogeneous/ww%g-kT%g.dat",
          lambda, temp);
  FILE *out = fopen(fname, "w");
  assert(out);
  delete[] fname;

  const double dff = 0.01;
  for (double ff = dff; ff < 0.7; ff += dff) {
    HomogeneousSW_liquid hf;
    hf.R() = radius;
    hf.epsilon() = epsilon;
    hf.kT() = temp;
    hf.lambda() = lambda;
    hf.n() = ff/(M_PI*uipow(sigma, 3)/6);
    hf.mu() = 0;
    hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf

    double pressure = hf.n()*hf.d_by_dn() - hf.energy();
    printf("%10g %10g %10g\n", ff, hf.energy(), pressure);
    fprintf(out, "%g\t%g\n", ff, pressure);
    // hf.printme("\t\t");
  }

  fclose(out);
  return 0;
}
