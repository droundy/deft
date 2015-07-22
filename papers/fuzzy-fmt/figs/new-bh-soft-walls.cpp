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
#include "new/WhiteBearFluidFast.h"
#include "new/WhiteBearFluidVeffFast.h"
#include "new/HomogeneousWhiteBearFluidFast.h"
#include "new/Minimize.h"

// Here we set up the lattice.
static double width = 30;
const double dx = 0.01;
const double dw = 0.01;
const double spacing = 3.0; // space on each side
const int N = 1000000;

double rad = 0;
const double rho = 1.0;
const double eps = 1.0;
const double sigma = 1.0; /* Let's define sigma == 1 for this one? */
const double sig6 = uipow(sigma,6);
const double sig12 = uipow(sig6,2);
const double R = sigma*pow(2.0,1.0/6.0);

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  double peak = peak_memory()/1024.0/1024;
  printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  fflush(stdout);
  last_time = t;
}

double R_BH(const double kT) {
  double bh_diameter = 0;
  const double dr = R/N;
  const double beta = 1.0/kT;
  for (double r_cur=dr/2; r_cur < R; r_cur += dr) {
    bh_diameter += (1 - exp(-beta*(4*eps*(uipow(sigma/r_cur,12) - uipow(sigma/r_cur,6)) + eps)))*dr;
  }
  return bh_diameter/2;
}

void run_soft_walls(double reduced_density, WhiteBearFluidVeff *f, double kT, double rad_bh) {
  Minimize min(f);
  min.set_relative_precision(0);
  min.set_maxiter(100);
  min.set_miniter(9);
  min.precondition(true);

  printf("========================================\n");
  printf("| Working on rho* = %4g and kT = %4g |\n", reduced_density, kT);
  printf("========================================\n");
  while (min.improve_energy(quiet)) {
    // f->run_finite_difference_test("WhiteBear");
  }
  took("Doing the minimization");
  min.print_info();

  char *fname = new char[5000];
  mkdir("papers/fuzzy-fmt/figs/new-data", 0777); // make sure the directory exists
  snprintf(fname, 5000, "papers/fuzzy-fmt/figs/new-data/bh-soft-wall-%04.2f-%04.2f.dat", reduced_density, kT);
  FILE *o = fopen(fname, "w");
  if (!o) {
    fprintf(stderr, "error creating file %s\n", fname);
    exit(1);
  }
  delete[] fname;
  const int Nz = f->Nz();
  Vector rz = f->get_rz();
  Vector n = f->get_n();
  printf("multiplying by sigma = %g\n", sigma);
  for (int i=0;i<Nz/2;i++) {
    fprintf(o, "%g\t%g\n", (rz[i] - spacing + R)/sigma, n[i]*uipow(sigma, 3));
  }
  fclose(o);
}

/*--

for kT in np.arange(0.1, 10.05, 0.1):
  for rho in np.arange(0.1, 2.05, 0.1):
    self.rule('%s %g %g' % (exe, rho, kT),
              [exe],
              ["papers/fuzzy-fmt/figs/new-data/bh-soft-wall-%04.2f-%04.2f.dat" % (rho, kT)])

--*/

int main(int argc, char **argv) {
  double reduced_density, temp;
  if (argc != 3) {
    printf("usage: %s reduced_density kT\n", argv[0]);
    return 1;
  }
  sscanf(argv[1], "%lg", &reduced_density);
  sscanf(argv[2], "%lg", &temp);

  HomogeneousWhiteBearFluid hf;

  rad = R_BH(temp);
  printf("rad is %g\n", rad);

  hf.R() = rad;
  hf.kT() = temp;
  hf.n() = reduced_density/uipow(sigma,3);
  printf("dividing by sigma = %g\n", sigma);
  printf("eta is %g\n", hf.n()*uipow(rad,3)*M_PI*4/3);
  hf.mu() = 0;
  hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf
  printf("bulk energy is %g\n", hf.energy());
  printf("cell energy should be %g\n", hf.energy()*dw*dw*width);

  WhiteBearFluidVeff f(dw, dw, width + spacing, dx);
  f.R() = hf.R();
  f.kT() = hf.kT();
  f.mu() = hf.mu();
  f.Vext() = 0;
  f.Veff() = -temp*log(hf.n());

  {
    const int Ntot = f.Nx()*f.Ny()*f.Nz();
    const Vector rz = f.get_rz();
    const double Vmax = 500*temp;
    for (int i=0; i<Ntot; i++) {
      if ( fabs(rz[i]) <= spacing ) {
        f.Vext()[i] = Vmax;
      } else if ( fabs(rz[i]) < spacing + rad) {
        double dist = fabs(rz[i]) - spacing;
        double dist3 = dist*dist*dist;
        double Vsw = 2*M_PI*rho*eps*((dist3-rad*rad*rad)/6
                                     + 2*sig12*(1/uipow(dist, 9) - 1/uipow(rad,9))/45
                                     + (rad - dist)*(rad*rad/2 + sig6/uipow(rad,4) 
                                                     - 2*sig12/5/uipow(rad,10))
                                     + sig6*(1/(rad*rad*rad) - 1/(dist3))/3);
        if (Vsw > Vmax) {
          f.Vext()[i] = Vmax; // this is "infinity" for us
        } else {
          f.Vext()[i] = Vsw; 
          printf("At distance %g, V_sw is %g\n", dist, Vsw);
        }
      } else {
        f.Vext()[i] = 0;
      }
    }
  }
  /* for (double rd = 0.50; rd < 0.63; rd += 0.01) {
    hf.n() = rd/uipow(sigma,3);
    f.Veff() = -temp*log(hf.n());
    printf("  %g\t%g\t%g\n", rd, f.energy(), hf.energy()*dw*dw*width);
    }*/

  printf("my energy is %g\n", f.energy());

  run_soft_walls(reduced_density, &f, temp, rad);
  return 0;
}
