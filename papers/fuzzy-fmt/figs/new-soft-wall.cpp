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
static double width = 20;
const double dx = 0.001;
const double dw = 0.002;
const double spacing = 1.5; // space on each side
const double epsilon = 1.0;
const double radius = 1.0;

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  double peak = peak_memory()/1024.0/1024;
  printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  fflush(stdout);
  last_time = t;
}

double soft_wall_potential(double z, double Vcutoff) {
  const double R_0 = 2*radius;
  const double rho = 1.0; // wall density
  z = fabs(z);

  if (z >= (spacing + R_0)) return 0;
  if ( z < spacing ) return Vcutoff;

  z = z - spacing; // now we set z=0 to be overlap with the wall (and infinite potential)

  const double sig = radius*pow(2,5.0/6.0);
  const double sig6 = pow(sig,6);
  const double sig12 = pow(sig,12);
  const double z3 = pow(z,3);
  const double z9 = pow(z,9);
  const double R3 = pow(R_0,3);
  const double R9 = pow(R_0,9);

  double potential = 2*M_PI*rho*epsilon*((z3-R3)/6
                                       + 2*sig12*(1/z9 - 1/R9)/45
                                       + (R_0 - z)*(R_0*R_0/2 + sig6/pow(R_0,4)
                                                    - 2*sig12/5/pow(R_0,10))
                                       + sig6*(1/R3-1/z3)/3);
  if (potential < Vcutoff) return potential;
  return Vcutoff;
}

void run_walls(double reduced_density, SFMTFluidVeff *f, double kT) {
  Minimize min(f);
  min.set_relative_precision(0);
  min.set_maxiter(100);
  min.set_miniter(9);
  min.precondition(true);
  if (reduced_density == 0.4 && kT == 0.01) min.set_known_true_energy(-2.41098243257584e-07);

  printf("========================================\n");
  printf("| Working on rho* = %4g and kT = %4g |\n", reduced_density, kT);
  printf("========================================\n");
  while (min.improve_energy(quiet)) {
    //f->run_finite_difference_test("SFMT");
  }
  took("Doing the minimization");
  min.print_info();

  char *fname = new char[5000];
  mkdir("papers/fuzzy-fmt/figs/new-data", 0777); // make sure the directory exists
  snprintf(fname, 5000, "papers/fuzzy-fmt/figs/new-data/soft-wall-%04.2f-%04.2f.dat", reduced_density, kT);
  FILE *o = fopen(fname, "w");
  if (!o) {
    fprintf(stderr, "error creating file %s\n", fname);
    exit(1);
  }
  delete[] fname;
  const int Nz = f->Nz();
  Vector Vext = f->Vext();
  Vector rz = f->get_rz();
  Vector n = f->get_n();
  Vector n3 = f->get_n3();
  for (int i=Nz/2+1;i<Nz;i++) {
    fprintf(o, "%g\t%g\t%g\n", rz[i] - spacing, n[i]/pow(2,-5.0/2.0), Vext[i]/kT);
  }
  for (int i=0;i<Nz/2;i++) {
    fprintf(o, "%g\t%g\t%g\n", rz[i] - spacing, n[i]/pow(2,-5.0/2.0), Vext[i]/kT);
  }
  fclose(o);
}

int main(int argc, char **argv) {
  double reduced_density, temp;
  if (argc != 3) {
    printf("usage: %s reduced_density kT\n", argv[0]);
    return 1;
  }
  sscanf(argv[1], "%lg", &reduced_density);
  sscanf(argv[2], "%lg", &temp);

  HomogeneousSFMTFluid hf;
  const double rad = 1;
  hf.sigma() = rad*pow(2,5.0/6.0);
  hf.epsilon() = 1;
  hf.kT() = temp;
  hf.n() = reduced_density*pow(2,-5.0/2.0);
  hf.mu() = 0;
  hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf
  printf("bulk energy is %g\n", hf.energy());
  //hf.printme("XXX:");
  printf("cell energy should be %g\n", hf.energy()*dw*dw*width);

  SFMTFluidVeff f(dw, dw, width, dx);
  f.sigma() = hf.sigma();
  f.epsilon() = hf.epsilon();
  f.kT() = hf.kT();
  //f.Veff() = 0;
  f.mu() = hf.mu();
  f.Vext() = 0;

  {
    const int Ntot = f.Nx()*f.Ny()*f.Nz();
    const Vector rz = f.get_rz();
    for (int i=0; i<Ntot; i++) {
      f.Vext()[i] = soft_wall_potential(rz[i], 20*temp);
      if (fabs(rz[i]) < spacing) {
        f.Veff()[i] = -temp*log(0.01*hf.n());
        //f.n()[i] = 0.01*hf.n();
      } else {
        f.Veff()[i] = -temp*log(hf.n());
        //f.n()[i] = hf.n();
      }
    }
  }
  printf("my energy is %g\n", f.energy());

  run_walls(reduced_density, &f, temp);
  return 0;
}
