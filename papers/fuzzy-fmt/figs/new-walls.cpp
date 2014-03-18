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
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include "new/SFMTFluidFast.h"
#include "new/SFMTFluidVeffFast.h"
#include "new/HomogeneousSFMTFluidFast.h"
#include "new/Minimize.h"

// Here we set up the lattice.
static double width = 15;
const double dx = 0.001;
const double dw = 0.002;
const double spacing = 1.5; // space on each side

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  double peak = peak_memory()/1024.0/1024;
  printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  fflush(stdout);
  last_time = t;
}

void run_walls(double eta, SFMTFluidVeff *f, double kT) {
  Minimize min(f);
  min.set_relative_precision(1e-12);
  min.set_maxiter(250);
  min.set_miniter(6);

  printf("======================================\n");
  printf("| Working on eta = %g and kT = %g |\n", eta, kT);
  printf("======================================\n");
  while (min.improve_energy(verbose)) {
    //f->run_finite_difference_test("SFMT");
  }
  took("Doing the minimization");
  min.print_info();

  char *fname = new char[5000];
  mkdir("papers/fuzzy-fmt/figs/new-data", 0777); // make sure the directory exists
  snprintf(fname, 5000, "papers/fuzzy-fmt/figs/new-data/wall-%04.2f-%08.5g.dat", eta, kT);
  FILE *o = fopen(fname, "w");
  if (!o) {
    fprintf(stderr, "error creating file %s\n", fname);
    exit(1);
  }
  delete[] fname;
  const int Nz = f->Nz();
  Vector rz = f->get_rz();
  Vector n = f->get_n();
  Vector n3 = f->get_n3();
  for (int i=0;i<Nz/2;i++) {
    fprintf(o, "%g\t%g\t%g\n", rz[i] - spacing, n[i]*(4*M_PI/3), n3[i]);
  }
  fclose(o);
}

int main(int, char **) {
  FILE *fout = fopen("papers/fuzzy-fmt/figs/wallsfillingfracInfo.txt", "w");
  fclose(fout);
  const double temps[] = { 0.01, 0.02, 0.03 };
  for (double eta = 0.4; eta > 0.05; eta-=0.1) {
    for (unsigned int i = 0; i<sizeof(temps)/sizeof(temps[0]); i++) {
      const double temp = temps[i];
      HomogeneousSFMTFluid hf;
      hf.R() = 1;
      hf.V0() = 1;
      hf.kT() = temp;
      hf.n() = eta/(4*M_PI/3);
      hf.mu() = 0;
      hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf
      printf("bulk energy is %g\n", hf.energy());
      //hf.printme("XXX:");
      printf("cell energy should be %g\n", hf.energy()*dw*dw*width);

      SFMTFluidVeff f(dw, dw, width + spacing, dx);
      f.R() = hf.R();
      f.V0() = hf.V0();
      f.kT() = hf.kT();
      f.Veff() = 0;
      f.mu() = hf.mu();
      f.Vext() = 0;

      {
        const int Ntot = f.Nx()*f.Ny()*f.Nz();
        const Vector rz = f.get_rz();
        for (int i=0; i<Ntot; i++) {
          if (fabs(rz[i]) < spacing) {
            f.Vext()[i] = 10*temp; // this is "infinity" for our wall
            f.Veff()[i] = -temp*log(hf.n());
          } else {
            f.Veff()[i] = -temp*log(hf.n());
          }
        }
      }
      printf("my energy is %g\n", f.energy());

      run_walls(eta, &f, temp);
    }
  }
  return 0;
}
