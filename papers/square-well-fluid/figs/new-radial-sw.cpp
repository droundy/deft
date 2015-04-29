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
#include "new/SW_liquidFast.h"
#include "new/HomogeneousSW_liquidFast.h"
#include "new/SW_liquidVeffFast.h"
#include "new/Minimize.h"
#include "version-identifier.h"

// Here we set up the lattice.
double zmax = 20;
double ymax = zmax;
double xmax = zmax;
double dx = 0.05;
const double epsilon = 1.0;
const double radius = 1.0;
const double sigma = 2*radius;
double lambda = 1;

static void took(const char *name) {
  static clock_t last_time = 0;
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  printf("\t\t%s took %g seconds\n", name, (t-last_time)/double(CLOCKS_PER_SEC));
  fflush(stdout);
  last_time = t;
}

void run_sw_liquid(double ff, SW_liquidVeff *f, double kT) {
  Minimize min(f);
  min.set_relative_precision(0);
  min.set_maxiter(100);
  min.precondition(true);

  char *dumpname = new char[5000];
  snprintf(dumpname, 5000, "papers/square-well-fluid/data/radial-sw-%.2f-%.2f-%.2f-X.dat",
           kT, f->lambda(), ff);
  f->get_rx().dumpSliceZ(dumpname, f->Nx(), f->Ny(), f->Nz(), 0);
  snprintf(dumpname, 5000, "papers/square-well-fluid/data/radial-sw-%.2f-%.2f-%.2f-Y.dat",
           kT, f->lambda(), ff);
  f->get_ry().dumpSliceZ(dumpname, f->Nx(), f->Ny(), f->Nz(), 0);
  snprintf(dumpname, 5000, "papers/square-well-fluid/data/radial-sw-%.2f-%.2f-%.2f-n.dat",
           kT, f->lambda(), ff);

  char *fname = new char[5000];
  //mkdir("papers/squre-well-fluid/figs/new-data", 0777); // make sure the directory exists
  snprintf(fname, 5000, "papers/square-well-fluid/data/radial-sw-%.2f-%.2f-%.2f.dat",
           kT, f->lambda(), ff);

  printf("========================================\n");
  printf("| Working on ff = %4g, lambda = %4g and kT = %4g |\n", ff, f->lambda(), kT);
  printf("========================================\n");
  while (min.improve_energy(verbose)) {

    took("Doing the minimization step");

    FILE *o = fopen(fname, "w");
    if (!o) {
      fprintf(stderr, "error creating file %s\n", fname);
      exit(1);
    }
    const int Nz = f->Nz();
    Vector Vext = f->Vext();
    Vector r = f->get_r();
    Vector n = f->get_n();
    for (int i=0;i<Nz/2;i++) {
      fprintf(o, "%g\t%g\t%g\n", r[i]/sigma, n[i]*uipow(sigma, 3)/6, Vext[i]);
    }
    fclose(o);

    n.dumpSliceZ(dumpname, f->Nx(), f->Ny(), f->Nz(), 0);

    took("Outputting to file");
  }
  min.print_info();
  delete[] fname;
  delete[] dumpname;
}

int main(int argc, char **argv) {
  double ff, temp, lambda;
  if (argc != 4) {
    printf("usage: %s ff lambda kT\n", argv[0]);
    return 1;
  }
  printf("git version: %s\n", version_identifier());

  sscanf(argv[1], "%lg", &ff);
  sscanf(argv[2], "%lg", &lambda);
  sscanf(argv[3], "%lg", &temp);

  printf("ff %g lam %g temp %g\n", ff, lambda, temp);

  HomogeneousSW_liquid hf;
  hf.R() = radius;
  hf.epsilon() = epsilon;
  hf.kT() = temp;
  hf.lambda() = lambda;
  hf.n() = ff*6/uipow(sigma, 3);
  hf.mu() = 0;
  printf("n3 = %g comes from n = %g\n", hf.get_n3(), hf.n());
  printf("bulk energy is not %g\n", hf.energy());
  hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf
  printf("bulk energy is %g\n", hf.energy());
  //hf.printme("XXX:");
  printf("cell energy should be %g\n", hf.energy()*xmax*ymax*zmax);

  SW_liquidVeff f(xmax, ymax, zmax, dx);
  f.R() = hf.R();
  f.epsilon() = hf.epsilon();
  f.lambda() = hf.lambda();
  f.kT() = hf.kT();
  //f.Veff() = 0;
  f.mu() = hf.mu();
  f.Vext() = 0;
  f.Veff() = -temp*log(hf.n()); // start with a uniform density as a guess

  {
    const int Ntot = f.Nx()*f.Ny()*f.Nz();
    const Vector r = f.get_r();
    for (int i=0; i<Ntot; i++) {
      const double Vmax = 100*temp;
      f.Vext()[i] = 0;
      if (r[i] > sigma && r[i] < sigma*lambda) {
        f.Vext()[i] = -epsilon;
      }
      if (r[i] <= sigma) {
        f.Vext()[i] = Vmax;
      }
    }
  }
  printf("my energy is %g\n", f.energy());
  took("Finding the energy a single time");

  run_sw_liquid(ff, &f, temp);
  return 0;
}
