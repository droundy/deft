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
#include "new/WhiteBearFluidFast.h"
#include "new/WhiteBearFluidVeffFast.h"
#include "new/HomogeneousWhiteBearFluidFast.h"
#include "new/Minimize.h"
#include "version-identifier.h"

// Here we set up the lattice.
double zmax = 16;
double ymax = zmax;
double xmax = zmax;
double dx = 0.05;
const double epsilon = 1.0;
const double radius = 1.0;
const double R = 2*radius;
const double sigma = R*pow(2,-1.0/6.0);
const int N = 1000000;

static void took(const char *name) {
  static clock_t last_time = 0;
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  printf("\t\t%s took %g seconds\n", name, (t-last_time)/double(CLOCKS_PER_SEC));
  fflush(stdout);
  last_time = t;
}

double R_BH(const double kT) {
  printf("kT for R_BH is %g.\n", kT);
  double bh_diameter = 0;
  const double dr = R/N;
  const double beta = 1.0/kT;
  printf("Beta is %g.\n", beta);
  for (double r_cur=dr/2; r_cur < R; r_cur += dr) {
    bh_diameter += (1 - exp(-beta*(4*epsilon*(uipow(sigma/r_cur,12)
                                              - uipow(sigma/r_cur,6)) + epsilon)))*dr;
  }
  return bh_diameter/2;
}

void run_minimization(double reduced_density,WhiteBearFluidVeff *f, double kT) {
  Minimize min(f);
  min.set_relative_precision(0);
  min.set_maxiter(1000);
  min.set_miniter(9);
  min.precondition(true);

  char *fname = new char[5000];
  mkdir("papers/fuzzy-fmt/figs/new-data", 0777); // make sure the directory exists
  snprintf(fname, 5000, "papers/fuzzy-fmt/figs/new-data/radial-bh-wca-%06.4f-%04.2f.dat", kT, reduced_density);

  printf("========================================\n");
  printf("| Working on rho* = %4g and kT = %4g |\n", reduced_density, kT);
  printf("========================================\n");
  do {
    //f->run_finite_difference_test("SFMT", 0, 100*min.recent_stepsize());

    took("Doing the minimization step");

    const int Nz = f->Nz();
    Vector Vext = f->Vext();
    Vector r = f->get_r();
    Vector n = f->get_n();
    f->get_Fideal(); // FIXME this is a hokey trick to make dV be defined
    Vector n3 = f->get_n3();

    FILE *o = fopen(fname, "w");
    if (!o) {
      fprintf(stderr, "error creating file %s\n", fname);
      exit(1);
    }
    for (int i=0;i<Nz/2;i++) {
      fprintf(o, "%g\t%g\t%g\t%g\n", r[i]/sigma, n[i]*uipow(sigma, 3), Vext[i], n3[i]);
    }
    fclose(o);

    took("Outputting to file");
  } while (min.improve_energy(gossipy));
  min.print_info();
  delete[] fname;
}

/* - -
# To run this automagically in fac, remove spaces from the above line

for kT in np.arange(0.1, 2.05, 0.1):
  for rho in np.arange(0.1, 2.05, 0.1):
    self.rule('%s %g %g' % (exe, rho, kT),
              [exe],
              ["papers/fuzzy-fmt/figs/new-data/radial-bh-%04.2f-%04.2f.dat" % (rho, kT)])

- - */

int main(int argc, char **argv) {
  double reduced_density, temp;
  if (argc != 3) {
    printf("usage: %s reduced_density kT\n", argv[0]);
    return 1;
  }
  printf("git version: %s\n", version_identifier());

  sscanf(argv[1], "%lg", &reduced_density);
  sscanf(argv[2], "%lg", &temp);

  HomogeneousWhiteBearFluid hf;
  printf("dx is %g\n", dx);

  double rad_bh = R_BH(temp);
  printf("rad is %g\n", rad_bh);

  hf.R() = rad_bh;
  hf.kT() = temp;
  hf.n() = reduced_density*pow(2,-5.0/2.0);
  printf("dividing by sigma = %g\n", sigma);
  printf("eta is %g\n", hf.n()*uipow(radius,3)*M_PI*4/3);
  hf.mu() = 0;
  hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf
  printf("bulk energy is %g\n", hf.energy());
  printf("cell energy should be %g\n", hf.energy()*dx*dx*dx);

  WhiteBearFluidVeff f(xmax, ymax, zmax, dx);
  f.R() = hf.R();
  f.kT() = hf.kT();
  f.mu() = hf.mu();
  f.Vext() = 0;
  f.Veff() = -temp*log(hf.n());

  {
    const int Ntot = f.Nx()*f.Ny()*f.Nz();
    const Vector r = f.get_r();
    for (int i=0; i<Ntot; i++) {
      const double Vmax = 100*temp;
      f.Vext()[i] = 4*epsilon*(uipow(sigma/r[i], 12) - uipow(sigma/r[i], 6)) + epsilon;
      if (r[i] > R) { f.Vext()[i] = 0; }
      if (!(f.Vext()[i] < Vmax)) f.Vext()[i] = Vmax;

      if (f.Vext()[i] > 0) {
        f.Veff()[i] += f.Vext()[i]; // adjust uniform guess based on repulsive potential
      }
    }
  }
  took("setting up the potential and Veff");

  run_minimization(reduced_density, &f, temp);

  return 0;
}
