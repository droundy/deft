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
#include <time.h>
#include "equation-of-state.h"
#include "LineMinimizer.h"
#include "ContactDensity.h"
#include "utilities.h"
#include "handymath.h"

#include "SoftFluidFast.h"
#include "HardRosenfeldFluidFast.h"

/*--

for kT in np.arange(0.1, 2.05, 0.1):
  for rho in np.arange(0.1, 2.05, 0.1):
    self.rule('%s %g %g' % (exe, rho, kT),
              [exe],
              ["papers/fuzzy-fmt/figs/wallssoft-%06.4f-%04.2f.dat" % (kT, rho)])

--*/

// Here we set up the lattice.
static double width = 15;
const double dx = 0.01;
const double dw = 0.01;
const double spacing = 1.5; // space on each side

double notinwall(Cartesian r) {
  const double z = r.z();
  if (fabs(z) > spacing) {
    return 1;
  }
  return 0;
}

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  double peak = peak_memory()/1024.0/1024;
  printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  fflush(stdout);
  last_time = t;
}

const int numiters = 25;


void z_plot(const char *fname, const Grid &a) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = 0;
  const int y = 0;
  for (int z=0; z<gd.Nz/2; z++) {
    Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
    double ahere = a(x,y,z);
    fprintf(out, "%g\t%g\n", here[2] - spacing, ahere);
  }
  fclose(out);
}

double run_walls(double reduced_density, const char *name, Functional fhs, double teff) {
  double kT = teff;
  if (kT == 0) kT = 1;

  Functional f = OfEffectivePotential(fhs);

  const double zmax = width + 2*spacing;
  Lattice lat(Cartesian(dw,0,0), Cartesian(0,dw,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, dx);

  Grid constraint(gd);
  constraint.Set(notinwall);
  f = constrain(constraint, f);

  Grid potential(gd);
  potential = pow(2,-5.0/2.0)*(reduced_density*constraint + 1e-4*reduced_density*VectorXd::Ones(gd.NxNyNz));
  potential = -kT*potential.cwise().log();

  const double approx_energy = fhs(kT, reduced_density*pow(2,-5.0/2.0))*dw*dw*width;
  const double precision = fabs(approx_energy*1e-11);
  printf("\tMinimizing to %g absolute precision from %g from %g...\n", precision, approx_energy, kT);
  fflush(stdout);

  Minimizer min = Precision(precision,
                            PreconditionedConjugateGradient(f, gd, kT,
                                &potential,
                                QuadraticLineMinimizer));
  took("Setting up the variables");
  if (strcmp(name, "hard") != 0 && false) {
    printf("For now, SoftFluid doesn't work properly, so we're skipping the\n");
    printf("minimization at temperature %g.\n", teff);
  } else {
    for (int i=0; min.improve_energy(false) && i<100; i++) {
    }
  }
  took("Doing the minimization");
  min.print_info();

  Grid density(gd, EffectivePotentialToDensity()(kT, gd, potential));
  //printf("# per area is %g at filling fraction %g\n", density.sum()*gd.dvolume/dw/dw, eta);

  char *plotname = (char *)malloc(1024);

  sprintf(plotname, "papers/fuzzy-fmt/figs/walls%s-%06.4f-%04.2f.dat", name, teff, reduced_density);
  z_plot(plotname, Grid(gd, density*pow(2,5.0/2.0)));
  free(plotname);

  took("Plotting stuff");
  printf("density %g gives ff %g for reduced density = %g and T = %g\n", density(0,0,gd.Nz/2),
         density(0,0,gd.Nz/2)*4*M_PI/3, reduced_density, teff);
  return density(0, 0, gd.Nz/2)*4*M_PI/3; // return bulk filling fraction
}

int main(int argc, char **argv) {
  double n_reduced, temp;
  if (argc != 3) {
    printf("usage: %s reduced_density kT\n", argv[0]);
    return 1;
  }
  sscanf(argv[1], "%lg", &n_reduced);
  sscanf(argv[2], "%lg", &temp);

  const double rad = 1; // //radius of our spheres
  const double sigma = rad*pow(2,5.0/6.0);

  Functional f = HardRosenfeldFluid(rad,0);
  if (temp > 0) f = SoftFluid(sigma, 1, 0);
  const double mu = find_chemical_potential(OfEffectivePotential(f), (temp)?temp:1, n_reduced*pow(2,-5.0/2.0));
  printf("mu is %g for reduced density = %g at temperature %g\n", mu, n_reduced, temp);
  if (temp > 0) f = SoftFluid(sigma, 1, mu);
  else f = HardRosenfeldFluid(rad, mu);

  const char *name = "hard";
  if (temp > 0) name = "soft";

  run_walls(n_reduced, name, f, temp);
  return 0;
}
