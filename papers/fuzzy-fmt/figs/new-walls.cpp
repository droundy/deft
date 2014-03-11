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

void run_walls(double eta, NewFunctional *f, double kT) {
  Minimize min(f);
  min.set_relative_precision(0);
  min.set_maxiter(25);

  while (min.improve_energy(min_details)) {
  }
  took("Doing the minimization");
  min.print_info();
  exit(1);
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
      // set mu based on derivative of hf
      const double derivative = hf.grad()[2];
      hf.mu() = derivative; // set eta to be the equilibrium density
      printf("bulk energy is %g\n", hf.energy());
      hf.printme("XXX:\t");
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
            f.Vext()[i] = 100*temp; // this is "infinity" for our wall
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
  // Just create this file so make knows we have run.
  if (!fopen("papers/fuzzy-fmt/figs/new-walls.dat", "w")) {
    printf("Error creating walls.dat!\n");
    return 1;
  }
  return 0;
}
