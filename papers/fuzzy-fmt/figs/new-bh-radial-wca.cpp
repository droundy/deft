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

#include "findbh.h"

// Here we set up the lattice.
double zmax = 16;
double ymax = zmax;
double xmax = zmax;
double dx = 0.02;

static void took(const char *name) {
  static clock_t last_time = 0;
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  printf("\t\t%s took %g seconds\n", name, (t-last_time)/double(CLOCKS_PER_SEC));
  fflush(stdout);
  last_time = t;
}

void run_minimization(double reduced_density,WhiteBearFluidVeff *f, double kT) {
  Minimize min(f);
  min.set_relative_precision(0);
  min.set_maxiter(1000);
  min.set_miniter(9);
  min.precondition(true);

  char *fname = new char[5000];
  mkdir("new-data", 0777); // make sure the directory exists
  snprintf(fname, 5000, "new-data/radial-bh-wca-%06.4f-%04.2f.dat", kT, reduced_density);

  printf("========================================\n");
  printf("| Working on rho* = %4g and kT = %4g |\n", reduced_density, kT);
  printf("========================================\n");
  do {
    //f->run_finite_difference_test("BH", 0, 100*min.recent_stepsize());

    took("Doing the minimization step");

    const int Nz = f->Nz();
    Vector r = f->get_r();
    Vector n = f->get_n();
    f->get_Fideal(); // FIXME this is a hokey trick to make dV be defined

    printf("The length of the vectors is %d", n.get_size());
    FILE *o = fopen(fname, "w");
    if (!o) {
      fprintf(stderr, "error creating file %s\n", fname);
      exit(1);
    }
    for (int i=0; i<Nz/2; i++) {
      fprintf(o, "%g\t%g\n", r[i]/sigma, n[i]*uipow(sigma, 3));
      if (i == int(Nz/2) - 1) {
        printf("\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n");
        printf("   n -> %.9g (i.e. g -> %g)\n",
           n[i]*uipow(sigma, 3), n[i]*uipow(sigma, 3)/reduced_density);
        printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
      }
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
              ["new-data/radial-bh-%04.2f-%04.2f.dat" % (rho, kT)])

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

  HomogeneousWhiteBearFluid hf = bh_homogeneous(reduced_density, temp);
  hf.mu() = 0;
  hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf
  printf("bulk energy is %g\n", hf.energy());
  //hf.printme("XXX:");
  printf("cell energy should be %g\n", hf.energy()*xmax*ymax*zmax);

  WhiteBearFluidVeff f = bh_inhomogeneous(temp, xmax, ymax, zmax, dx);  
  f.mu() = hf.mu();
  f.Vext() = 0;
  f.Veff() = -temp*log(hf.n()); // start with a uniform density as a guess
  {
    const int Ntot = f.Nx()*f.Ny()*f.Nz();
    const Vector r = f.get_r();
    for (int i=0; i<Ntot; i++) {
      const double Vmax = 100*temp;
      f.Vext()[i] = Vwca(r[i]);
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
