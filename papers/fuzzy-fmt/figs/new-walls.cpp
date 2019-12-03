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

#include "findxi.h"

// Here we set up the lattice.
static double width = 30;
const double dx = 0.01;
const double dw = 0.01;
const double spacing = 3.0; // space on each side

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  double peak = peak_memory()/1024.0/1024;
  printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  fflush(stdout);
  last_time = t;
}

void run_walls(double reduced_density, SFMTFluidVeff *f, double kT) {
  Minimize min(f);
  min.set_relative_precision(0);
  min.set_maxiter(10000);
  min.set_miniter(9);
  min.precondition(true);
  if (reduced_density == 0.4 && kT == 0.01) min.set_known_true_energy(-2.41098243257584e-07);

  {
    Vector phi3 = f->get_phi3_density();
    Vector rz = f->get_rz();
    Vector n3 = f->get_n3();
    Vector n = f->get_n();
    Vector trace_tensor_cubed = f->get_trace_tensor_cubed();
    Vector tensor_vector = f->get_tensor_vector();
    Vector n2x = f->get_n2vx();
    Vector n2y = f->get_n2vy();
    Vector n2z = f->get_n2vz();
    Vector n2xx = f->get_n2xx();
    Vector n2yy = f->get_n2yy();
    Vector n2zz = f->get_n2zz();
    Vector n2xy = f->get_n2xy();
    Vector n2yz = f->get_n2yz();
    Vector n2zx = f->get_n2zx();
    for (int i=0; i<phi3.get_size(); i++) {
      printf("   %8g %15g %15g %15g %15g %15g\n", rz[i], n[i], trace_tensor_cubed[i], n2xx[i], n2yy[i], n2zz[i]);
      if (rz[i] == rz[0] && i > 0) break;
    }
    exit(1);
  }
  printf("========================================\n");
  printf("| Working on rho* = %4g and kT = %4g |\n", reduced_density, kT);
  printf("========================================\n");
  while (min.improve_energy(verbose)) {
    // f->run_finite_difference_test("SFMT");
    Vector phi3 = f->get_phi3_density();
    Vector rz = f->get_rz();
    Vector n3 = f->get_n3();
    Vector n = f->get_n();
    Vector trace_tensor_cubed = f->get_trace_tensor_cubed();
    Vector tensor_vector = f->get_tensor_vector();
    for (int i=0; i<phi3.get_size(); i++) {
      printf("   %8g %15g %15g %15g %15g\n", rz[i], n[i], phi3[i], trace_tensor_cubed[i], tensor_vector[i]);
      if (rz[i] == rz[0] && i > 0) break;
    }
    if (min.energy() < -1e20) {
      printf("trouble with the energy: %g\n", min.energy());
      break;
    }
  }
  took("Doing the minimization");
  min.print_info();

  char *fname = new char[5000];
  mkdir("papers/fuzzy-fmt/figs/new-data", 0777); // make sure the directory exists
  snprintf(fname, 5000, "papers/fuzzy-fmt/figs/new-data/wall-%04.2f-%04.2f.dat", reduced_density, kT);
  FILE *o = fopen(fname, "w");
  if (!o) {
    fprintf(stderr, "error creating file %s\n", fname);
    exit(1);
  }
  delete[] fname;
  const int Nz = f->Nz();
  Vector rz = f->get_rz();
  Vector n = f->get_n();
  for (int i=0; i<Nz/2; i++) {
    fprintf(o, "%g\t%g\n", (rz[i] - spacing)/f->sigma(), n[i]*uipow(f->sigma(), 3));
  }
  fclose(o);
}

/*--

for kT in np.arange(0.1, 10.05, 0.1):
  for rho in np.arange(0.1, 2.05, 0.1):
    self.rule('%s %g %g' % (exe, rho, kT),
              [exe],
              ["papers/fuzzy-fmt/figs/new-data/wall-%04.2f-%04.2f.dat" % (rho, kT)])

--*/

int main(int argc, char **argv) {
  double reduced_density, temp;
  if (argc != 3) {
    printf("usage   : %s reduced_density kT\n", argv[0]);
    return 1;
  }
  sscanf(argv[1], "%lg", &reduced_density);
  sscanf(argv[2], "%lg", &temp);

  HomogeneousSFMTFluid hf = sfmt_homogeneous(reduced_density, temp);   //note: homogeneousFE/atom does not depend on fv or gw
  hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf
  printf("bulk energy is %g\n", hf.energy());
  //hf.printme("XXX:");
  printf("cell energy should be %g\n", hf.energy()*dw*dw*width);
  hf.printme("homogeneous");

  SFMTFluidVeff f = sfmt_inhomogeneous(temp, dw, dw, width + spacing, dx);
  f.mu() = hf.mu();
  f.Vext() = 0;
  //f.n() = hf.n();
  f.Veff() = -temp*log(hf.n());

  {
    const int Ntot = f.Nx()*f.Ny()*f.Nz();
    const Vector rz = f.get_rz();
    for (int i=0; i<Ntot; i++) {
      if (fabs(rz[i]) < spacing) {
        const double Vmax = 500*temp;
        f.Vext()[i] = Vmax; // this is "infinity" for our wall
        f.Veff()[i] = Vmax;
      } else {
        f.Vext()[i] = 0;
      }
    }
  }
  printf("my energy is %g\n", f.energy());

  run_walls(reduced_density, &f, temp);
  return 0;
}
