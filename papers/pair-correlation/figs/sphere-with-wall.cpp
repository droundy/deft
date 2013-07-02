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
#include <math.h>
#include "OptimizedFunctionals.h"
#include "equation-of-state.h"
#include "LineMinimizer.h"
#include "ContactDensity.h"
#include "utilities.h"
#include "handymath.h"
#include "errno.h"
#include "sys/stat.h" // for mkdir

int count =0;

// Maximum and spacing values for plotting
const double zmax = 20;
const double xmax = 10;
const double dx = 0.1;


const int num_eta = 10;
const double eta_step = 0.5/num_eta;

// Here we set up the lattice.
static double width = 20;
const double dw = 0.001;
const double spacing = 3; // space on each side


double notinwall_or_sphere(Cartesian r) {
  const double x = r.x();
  const double y = r.y();
  const double z = r.z();
  if (fabs(z) < spacing) {
    return 0;
  }
  const double dis = sqrt((z-spacing)*(z-spacing)+y*y*+x*x);
  if (z>0 && dis<2) {
    return 0;
  }
  return 1;
}

void z_plot(const char *fname, const Grid &a) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = 0;
  for (int y=0; y<gd.Nz/2; y++) {
    for (int z=0; z<gd.Nz/2; z++) {
      Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
      double density = a(x,y,z);
      double pair_correlation = a(x,y,z)/a(x,y,-z);/*this is density
      on the sphere side divided by the density the same distance from
      the wall on the other side (where there is no sphere)*/
      fprintf(out, "%g\t%g\t%g\t%g\n", here[2], here[1], density,
              pair_correlation); }
  }
  fclose(out);
}

Functional WB = HardSpheresNoTensor2(1.0);

int main(int, char **) {
  for (double eta = 0.3; eta < 0.35; eta += 0.1) {
    // Generates a data file for the pair distribution function, for filling fraction eta
    // and distance of first sphere from wall of z0. Data saved in a table such that the
    // columns are x values and rows are z1 values.
    printf("Now starting sphere_with_wall with eta = %g\n",eta);

    Functional f = OfEffectivePotential(WB + IdealGas());
    double mu = find_chemical_potential(f, 1, eta/(4*M_PI/3));
    f = OfEffectivePotential(WB + IdealGas()
                             + ChemicalPotential(mu));

    Lattice lat(Cartesian(width,0,0), Cartesian(0,width,0), Cartesian(0,0,width+2*spacing));
    GridDescription gd(lat, 0.1); // the resolution here dramatically affects our memory use

    Grid potential(gd);
    Grid constraint(gd);

    constraint.Set(*notinwall_or_sphere);
    f = constrain(constraint, f);

    potential = (eta*constraint + 1e-4*eta*VectorXd::Ones(gd.NxNyNz))/(4*M_PI/3);
    potential = -potential.cwise().log();

    const double approx_energy = (WB + IdealGas() + ChemicalPotential(mu))(1, eta/(4*M_PI/3))*dw*dw*width;
    const double precision = fabs(approx_energy*1e-4);
    //printf("Minimizing to %g absolute precision...\n", precision);
    Minimizer min = Precision(precision,
                              PreconditionedConjugateGradient(f, gd, 1,
                                                              &potential,
                                                              QuadraticLineMinimizer));
    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      fflush(stdout);
    }

    for (int i=0;min.improve_energy(true) && i<100;i++) {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      fflush(stdout);
    }

    Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));

    char *plotname = new char[1024];
    sprintf(plotname, "papers/pair-correlation/figs/wallsWB-with-sphere-%04.2f.dat", eta);
    z_plot(plotname, density);
    delete[] plotname;

    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      fflush(stdout);
    }
  }

  // Just create this file so make knows we have run.
  if (!fopen("papers/pair-correlation/figs/walls_sphere.dat", "w")) {
    printf("Error creating walls.dat!\n");
    return 1;
  }
  return 1;
}
