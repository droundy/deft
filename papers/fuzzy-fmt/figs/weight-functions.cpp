// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2015 The Deft Authors
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
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include "new/HomogeneousSFMTFluidFast.h"
#include "new/SFMTFluidVeffFast.h"
#include "new/Minimize.h"
#include "version-identifier.h"

// Here we set up the lattice.
double zmax = 8;
double ymax = zmax;
double xmax = zmax;
double dx = 0.025;

#include "findxi.h"

static void took(const char *name) {
  static clock_t last_time = 0;
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  printf("\t\t%s took %g seconds\n", name, (t-last_time)/double(CLOCKS_PER_SEC));
  fflush(stdout);
  last_time = t;
}

int main(int argc, char **argv) {
  double temp;
  if (argc != 2) {
    printf("usage: %s kT\n", argv[0]);
    return 1;
  }
  printf("git version: %s\n", version_identifier());

  sscanf(argv[1], "%lg", &temp);

  SFMTFluidVeff f = sfmt_inhomogeneous(temp, xmax, ymax, zmax, dx);
  f.Vext() = 0;
  f.Veff() = -temp*log(0.00001); // essentially zero everywhere else
  f.Veff()[0] = temp*log(1.0/dx/dx/dx); // here is the delta function.

  took("initializing functional");
  const int Nz = f.Nz();
  {
    Vector n3 = f.get_n3();
    took("finding n3");
    Vector n2 = f.get_n2();
    took("finding n2");
    Vector n1 = f.get_n1();
    took("finding n1");
    Vector n0 = f.get_n0();
    took("finding n0");
    Vector rx = f.get_rx();
    Vector ry = f.get_ry();
    Vector rz = f.get_rz();
    Vector n = f.get_n();
    took("finding r vector and n");
    Vector n2x = f.get_n2vx();
    took("finding n2x");
    Vector n2y = f.get_n2vy();
    took("finding n2y");
    Vector n2z = f.get_n2vz();
    took("finding n2z");
    Vector n2xx = f.get_n2xx();
    took("finding n2xx");
    Vector n2yy = f.get_n2yy();
    took("finding n2yy");
    Vector n2zz = f.get_n2zz();
    took("finding n2zz");
    Vector n2xy = f.get_n2xy();
    took("finding n2xy");
    Vector n2yz = f.get_n2yz();
    took("finding n2yz");
    Vector n2zx = f.get_n2zx();
    took("finding n2zx");
    char *fname = malloc(4096);
    sprintf(fname, "weight-functions-%g.dat", temp);
    printf("creating file %s\n", fname);
    FILE *o = fopen(fname, "w");
    if (!o) {
      fprintf(stderr, "error creating file\n");
      exit(1);
    }
    fprintf(o, "# %s\t%s\t%s", "x", "y", "z");
    fprintf(o, "\t%s\t%s\t%s\t%s\t%s", "n", "n3", "n2", "n1", "n0");
    fprintf(o, "\t%s\t%s\t%s", "n2x", "n2y", "n2z");
    fprintf(o, "\t%s\t%s\t%s",  "n2xx", "n2yy", "n2zz");
    fprintf(o, "\t%s\t%s\t%s\n",  "n2xy", "n2yz", "n2zx");
    for (int i=0; i<Nz/2; i++) {
      fprintf(o, "%g\t%g\t%g", rx[i], ry[i], rz[i]);
      fprintf(o, "\t%g\t%g\t%g\t%g\t%g", n[i], n3[i], n2[i], n1[i], n0[i]);
      fprintf(o, "\t%g\t%g\t%g", n2x[i], n2y[i], n2z[i]);
      fprintf(o, "\t%g\t%g\t%g",  n2xx[i], n2yy[i], n2zz[i]);
      fprintf(o, "\t%g\t%g\t%g\n",  n2xy[i], n2yz[i], n2zx[i]);
    }
    fclose(o);
  }
  return 0;
}
