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
#include "version-identifier.h"

const bool use_veff = true;

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  double peak = peak_memory()/1024.0/1024;
  printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  fflush(stdout);
  last_time = t;
}

double inhomogeneity(Vector n) {
  double maxn = n[0];
  double minn = n[0];
  for (int i=0;i<n.get_size();i++) {
    if (n[i] > maxn) maxn = n[i];
    if (n[i] < minn) minn = n[i];
  }
  return (maxn - minn)/fabs(minn);
}

void run_solid(double lattice_constant, double reduced_density, double kT,
               SFMTFluidVeff *fveff, SFMTFluid *f,
               double homogeneous_free_energy) {
  Minimize min = (use_veff) ? Minimize(fveff) : Minimize(f);
  min.set_relative_precision(1e-12);
  min.set_maxiter(10000);
  min.set_miniter(9);
  min.precondition(true);

  printf("========================================\n");
  printf("| Working on rho* = %4g and kT = %4g and a = %g |\n", reduced_density, kT, lattice_constant);
  printf("========================================\n");
  while (min.improve_energy(verbose)) {
    //f->run_finite_difference_test("SFMT");
    printf("Compare with homogeneous free energy: %.15g\n", homogeneous_free_energy);
    Vector n = (use_veff) ? fveff->get_n() : f->n();
    double inh = inhomogeneity(n);
    printf("Inhomogeneity is %g\n", inh);
    double Ntot = n.sum()*f->get_dV();
    printf("Ntot = %g\n", Ntot);
    printf("n*bar = %g\n", Ntot/f->get_volume());
    printf("\n");
    if (inh < 1) {
      printf("It is flat enough for me!\n");
      break;
    }
  }
  took("Doing the minimization");
  min.print_info();

  char *fname = new char[5000];
  mkdir("papers/fuzzy-fmt/figs/new-data", 0777); // make sure the directory exists
  snprintf(fname, 5000, "papers/fuzzy-fmt/figs/new-data/melting-%04.2f-%04.2f-%04.2f.dat",
           lattice_constant, reduced_density, kT);
  FILE *o = fopen(fname, "w");
  if (!o) {
    fprintf(stderr, "error creating file %s\n", fname);
    exit(1);
  }
  delete[] fname;
  const int Nz = f->Nz();
  Vector rz = f->get_rz();
  Vector n = (use_veff) ? fveff->get_n() : f->n();
  for (int i=0;i<Nz/2;i++) {
    fprintf(o, "%g\t%g\n", rz[i], n[i]);
  }
  fclose(o);
}

int main(int argc, char **argv) {
  double lattice_constant, reduced_density, temp;
  if (argc != 4) {
    printf("usage: %s reduced-lattice-constant nliquid-reduced kT\n", argv[0]);
    return 1;
  }
  printf("git version: %s\n", version_identifier());
  sscanf(argv[1], "%lg", &lattice_constant);
  sscanf(argv[2], "%lg", &reduced_density);
  sscanf(argv[3], "%lg", &temp);

  HomogeneousSFMTFluid hf;
  hf.sigma() = 1;
  hf.epsilon() = 1;
  hf.kT() = temp;
  hf.n() = reduced_density;
  hf.mu() = 0;
  hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf

  const double homogeneous_free_energy = hf.energy()*lattice_constant*lattice_constant*lattice_constant;
  printf("bulk energy is %g\n", hf.energy());
  printf("liquid cell free energy should be %g\n", homogeneous_free_energy);

  const double dx = 0.05;
  SFMTFluidVeff fveff(lattice_constant, lattice_constant, lattice_constant, dx);
  SFMTFluid f(lattice_constant, lattice_constant, lattice_constant, dx);
  f.sigma() = hf.sigma();
  f.epsilon() = hf.epsilon();
  f.kT() = hf.kT();
  f.mu() = hf.mu();
  f.Vext() = 0;
  f.n() = hf.n();

  fveff.sigma() = hf.sigma();
  fveff.epsilon() = hf.epsilon();
  fveff.kT() = hf.kT();
  fveff.mu() = hf.mu();
  fveff.Vext() = 0;
  fveff.Veff() = 0;

  {
    // This is where we set up the inhomogeneous n(r)
    const int Ntot = f.Nx()*f.Ny()*f.Nz();
    const Vector rrx = f.get_rx();
    const Vector rry = f.get_ry();
    const Vector rrz = f.get_rz();
    const double gwidth = 0.1;
    const double norm = 1.2*pow(sqrt(2*M_PI)*gwidth, 3); // the prefactor is a safety factor
    Vector setn = (use_veff) ? fveff.Veff() : f.n();
    for (int i=0; i<Ntot; i++) {
      const double rx = rrx[i];
      const double ry = rry[i];
      const double rz = rrz[i];
      setn[i] = 0.0000001*hf.n();
      {
        double dist = sqrt(rx*rx + ry*ry+rz*rz);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;
      }
      {
        double dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                           (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                           rz*rz);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                    (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                    rz*rz);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                    (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                    rz*rz);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                    (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                    rz*rz);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;
      }
      {
        double dist = sqrt((rz-lattice_constant/2)*(rz-lattice_constant/2) +
                           (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                           rx*rx);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rz+lattice_constant/2)*(rz+lattice_constant/2) +
                    (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                    rx*rx);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rz-lattice_constant/2)*(rz-lattice_constant/2) +
                    (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                    rx*rx);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rz+lattice_constant/2)*(rz+lattice_constant/2) +
                    (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                    rx*rx);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;
      }
      {
        double dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                           (rz-lattice_constant/2)*(rz-lattice_constant/2) +
                           ry*ry);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                    (rz-lattice_constant/2)*(rz-lattice_constant/2) +
                    ry*ry);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                    (rz+lattice_constant/2)*(rz+lattice_constant/2) +
                    ry*ry);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;

        dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                    (rz+lattice_constant/2)*(rz+lattice_constant/2) +
                    ry*ry);
        setn[i] += exp(-0.5*dist*dist/gwidth/gwidth)/norm;
      }
    }
    if (use_veff) {
      for (int i=0; i<Ntot; i++) {
        fveff.Veff()[i] = -temp*log(fveff.Veff()[i]); // convert from density to effective potential
      }
    }
  }

  {
    char *fname = new char[5000];
    mkdir("papers/fuzzy-fmt/figs/new-data", 0777); // make sure the directory exists
    snprintf(fname, 5000, "papers/fuzzy-fmt/figs/new-data/initial-melting-%04.2f-%04.2f-%04.2f.dat",
             lattice_constant, reduced_density, temp);
    FILE *o = fopen(fname, "w");
    if (!o) {
      fprintf(stderr, "error creating file %s\n", fname);
      exit(1);
    }
    delete[] fname;
    const int Nz = f.Nz();
    Vector rz = f.get_rz();
    Vector n = (use_veff) ? fveff.get_n() : f.n();
    for (int i=0;i<Nz/2;i++) {
      fprintf(o, "%g\t%g\n", rz[i], n[i]);
    }
    fclose(o);
  }

  printf("my initial energy is %g\n", f.energy());
  if (f.energy() != f.energy()) {
    printf("FAIL!  nan for initial energy is bad!\n");
    exit(1);
  }

  run_solid(lattice_constant, reduced_density, temp, &fveff, &f, homogeneous_free_energy);
  return 0;
}
