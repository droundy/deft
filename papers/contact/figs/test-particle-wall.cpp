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
//#inc// lude <time.h>
#include "OptimizedFunctionals.h"
#include "equation-of-state.h"
#include "LineMinimizer.h"
#include "ContactDensity.h"
#include "utilities.h"
#include "handymath.h"

// Here we set up the lattice.
static const double width = 20;
const double dw = 0.001;
const double spacing = 3; // space on each sidestatic double diameter;
static double z_part;
static double diameter;

double notinwall(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  const double x = r.x();
  if (fabs(z) > spacing && sqrt(sqr(z-z_part)+sqr(y)+sqr(x)) > diameter) {
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
  last_time = t;
}

Functional WB = HardSpheresNoTensor(1.0);
Functional WBT = HardSpheresWBFast(1.0);
Functional WBm2 = HardSpheresWBm2(1.0);

const int numiters = 25;


void radial_plot(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d, const Grid &e,
                 const Grid &f, const Grid &g, const Grid &h) {
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
    double ahere = a(here);
    double bhere = b(here);
    double chere = c(here);
    double dhere = d(here);
    double ehere = e(here);
    double fhere = f(here);
    double ghere = g(here);
    double hhere = h(here);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[2],
            ahere, bhere, chere, dhere, ehere, fhere, ghere, hhere);
  }
  fclose(out);
}


const double golden = (1 + sqrt(5))/2;
const Cartesian icosohedron[] = {
  Cartesian( 1, 1, 1),
  Cartesian(-1, 1, 1),
  Cartesian( 1,-1, 1),
  Cartesian( 1, 1,-1),
  Cartesian( 1,-1,-1),
  Cartesian(-1, 1,-1),
  Cartesian(-1,-1, 1),
  Cartesian(-1,-1,-1),
  Cartesian(0, 1/golden, golden),
  Cartesian(0, 1/golden,-golden),
  Cartesian(0,-1/golden, golden),
  Cartesian(0,-1/golden,-golden),
  Cartesian( golden,0, 1/golden),
  Cartesian(-golden,0, 1/golden),
  Cartesian( golden,0,-1/golden),
  Cartesian(-golden,0,-1/golden),
  Cartesian( 1/golden, golden,0),
  Cartesian( 1/golden,-golden,0),
  Cartesian(-1/golden, golden,0),
  Cartesian(-1/golden,-golden,0),
};

// The following computes the closest value that is at least a vector
// delta from pos.
double closest_value(const Grid &g, Cartesian pos, Vector3d delta) {
  double closest = 1e300;
  double bestyet = 0;
  //printf("finding closest_value to %g %g %g\n", pos[0], pos[1], pos[2]);
  //printf("in direction %g %g %g\n", delta[0], delta[1], delta[2]);
  const GridDescription gd = g.description();
  for (int x=-gd.Nx/2; x<gd.Nx/2; x++) {
    for (int y=-gd.Ny/2; y<gd.Ny/2; y++) {
      for (int z=-gd.Nz/2; z<gd.Nz/2; z++) {
        Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
        //printf("Hello %d %d %d vs %d %d %d\n", x, y, z, gd.Nx, gd.Ny, gd.Nz);
        //printf("here is %g %g %g\n", here[0], here[1], here[2]);
        if ((here-pos-delta).dot(delta) > 0 && (here-pos).norm() < closest) {
          //printf("here we are... %d %d %d\n", x, y, z);
          bestyet = g(here);
          closest = (here-pos).norm();
        }
        //printf("Got here now..;\n");
      }
    }
  }
  return bestyet;
}
// The following computes the average value for the grid x at radius
// from pos.
void radial_average(const char *fname, const Grid &g, Cartesian pos, double radius) {
  FILE *out = fopen(fname, "w");
  if (!out){
    fprintf(stderr, "Unable to create file %s!\n", fname);
    return;
  }
  double total = 0;
  for (int i=0;i<20;i++) {
    total += closest_value(g, pos, Vector3d(radius * icosohedron[i]/icosohedron[i].norm()));
  }
  double z_part = pos[2];
  fprintf(out, "%g\t%g\n", z_part, total/20);
  return;
}

void plane_plot(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d, const Grid &e,
                 const Grid &f, const Grid &g, const Grid &h) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = 0;
  for (int z=0; z<gd.Nz/2; z++) {
    for (int y=-gd.Ny/2; y<gd.Ny/2; y++){
      Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
      double ahere = a(here);
      double bhere = b(here);
      double chere = c(here);
      double dhere = d(here);
      double ehere = e(here);
      double fhere = f(here);
      double ghere = g(here);
      double hhere = h(here);
      fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[2], here[1],
              ahere, bhere, chere, dhere, ehere, fhere, ghere, hhere);
    }
  }
  fclose(out);
}

void run_spherical_solute(double diam, double eta, double z_particle, const char *name, Functional fhs) {
  diameter = diam;
  z_part = z_particle;
  printf("Diameter is %g hard spheres, filling fraction %g, z position of particle %g\n", diameter, eta, z_part);

  const double meandensity = eta/(4*M_PI/3);

  Functional f = OfEffectivePotential(HardSpheresNoTensor(1.0) + IdealGas());
  double mu = find_chemical_potential(f, 1, meandensity);
  f = OfEffectivePotential(HardSpheresNoTensor(1.0) + IdealGas()
                           + ChemicalPotential(mu));

  const double zmax = width + 2*spacing;
  const double xmax = width;
  Lattice lat(Cartesian(xmax,0,0), Cartesian(0,xmax,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, 0.1);

  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
  took("Setting the constraint");

  f = constrain(constraint, f);

  constraint.epsNativeSlice("constraint-NativeSlice.eps", Cartesian(0,2*xmax,0), Cartesian(0,0,zmax), Cartesian(0,-xmax,0));

  potential = meandensity*constraint + 1e-4*meandensity*VectorXd::Ones(gd.NxNyNz);
  potential = -potential.cwise().log();
  //f.run_finite_difference_test("foobar", 1, potential);

  {
    const double approx_energy = (fhs + IdealGas() + ChemicalPotential(mu))(1, eta/(4*M_PI/3))*xmax*xmax*xmax;
    const double precision = fabs(approx_energy*1e-10);
    printf("Minimizing to %g absolute precision...\n", precision);
    Minimizer min = Precision(precision,
                              PreconditionedConjugateGradient(f, gd, 1,
                                                              &potential,
                                                              QuadraticLineMinimizer));
    for (int i=0;min.improve_energy(true) && i<200;i++) {
    }

    double energy = min.energy();
    printf("Energy is %.15g\n", energy);
    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
    }

    Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));

    printf("%g\t%.15g\n", diameter, energy);
  }

  Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }
  char *contact_plot = (char *)malloc(1024);
  sprintf(contact_plot, "papers/contact/figs/test-particle-cont-dens-%04.2f-%04.2f.dat",
          z_part, eta);
  radial_average(contact_plot, density, Cartesian(0,0,z_part), 2.0);
//printf("Density at contact is %g\n", radial_average(density, Cartesian(0,0,z_part), 1.0));
  printf("N = %g\n", density.sum()*gd.dvolume);
  char *plotname = (char *)malloc(1024);
  char *plane_plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/contact/figs/test-particle-wall%s-%04.1f-%04.2f-%04.2f.dat", name, diameter, eta, z_part);
  printf("Saving as %s\n", plotname);
  sprintf(plane_plotname, "papers/contact/figs/test-particle-wall%s-twodim-%04.1f-%04.2f-%04.2f.dat", name, diameter, eta, z_part);
  printf("Saving as %s\n", plane_plotname);
  if (plane_plotname == 0){
    printf("plane_plotname has not been allocated correctly!\n");
    exit (1);
  }
  Grid energy_density(gd, f(1, gd, potential));
  Grid correlation_S(gd, Correlation_S(1.0)(1, gd, density));
  Grid correlation_A(gd, Correlation_A(1.0)(1, gd, density));
  if (strlen(name) == 4) {
    correlation_S = Correlation_S_WBm2(1.0)(1, gd, density);
    correlation_A = Correlation_A_WBm2(1.0)(1, gd, density);
  }
  Grid gross_correlation(gd, GrossCorrelation(1.0)(1, gd, density));
  Grid n0(gd, ShellConvolve(1)(1, density)/(4*M_PI));
  Grid nA(gd, ShellConvolve(2)(1, density)/(4*M_PI*4));
  Grid yuwu_correlation(gd, YuWuCorrelation_S(1.0)(1, gd, density));
  sprintf(plotname, "papers/contact/figs/test-particle-wall%s-%04.1f-%04.2f-%04.2f.dat", name, diameter, eta, z_part);
  sprintf(plane_plotname, "papers/contact/figs/test-particle-wall%s-twodim-%04.1f-%04.2f-%04.2f.dat", name, diameter, eta, z_part);
  radial_plot(plotname, density, energy_density, correlation_S, yuwu_correlation,
              correlation_A, n0, gross_correlation, nA);
  plane_plot(plane_plotname, density, energy_density, correlation_S, yuwu_correlation,
              correlation_A, n0, gross_correlation, nA);
  free(plotname);
  free(plane_plotname);
  //  free(plane_plotname);

  {
    const GridDescription gdp = density.description();
    double inner_rad = diameter/2.0;

    double Ntot = density.sum()*gdp.dvolume;
    double Ndisplaced = eta*gdp.Lat.volume()/(4*M_PI/3) - Ntot;

    double mc_side_len = 25;
    double N = eta*mc_side_len*mc_side_len*mc_side_len/(4.0/3.0*M_PI) - Ndisplaced;

    FILE *fout = fopen("papers/contact/figs/testfillingfracInfo.txt", "a");
    if (fout==0){
      printf("Not able to open papers/contact/figs/testfillingfracInfo.txt properly\n");
      exit (1);
    }

    fprintf(fout, "For filling fraction %04.02f, inner-sphere size %04.02f and walls of length %04.02f you'll want to use %.0f spheres.\n\n", eta, inner_rad, mc_side_len, N);

    fclose(fout);
  }



  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }

  took("Plotting stuff");
}

int main(int argc, char *argv[]) {
  double eta;
  if (argc != 5) {
    printf("got argc %d\n", argc);
    printf("usage: %s diameter filling_fraction z_part (WB|WBT|WBm2)\n", argv[0]);
    return 1;
  }
  if (sscanf(argv[1], "%lg", &diameter) != 1) {
    printf("Got bad diameter argument: %s\n", argv[1]);
    return 1;
  }
  if (sscanf(argv[2], "%lg", &eta) != 1) {
    printf("Got bad eta argument: %s\n", argv[2]);
    return 1;
  }
  if (sscanf(argv[3], "%lg", &z_part) != 1) {
    printf("Got bad z_part argument: %s\n", argv[3]);
    return 1;
  }

  FILE *fout = fopen("papers/contact/figs/testfillingfracInfo.txt", "w");
  if (fout==0){
    printf("Not able to open papers/contact/figs/testfillingfracInfo.txt properly\n");
    exit (1);
  }
  fclose(fout);
  if (strlen(argv[4]) == 2) {
    run_spherical_solute(diameter, eta, z_part, "WB", WB);
  } else if (strlen(argv[4]) == 3) {
    run_spherical_solute(diameter, eta, z_part, "WBT", WBT);
  } else if (strlen(argv[4]) == 4) {
    run_spherical_solute(diameter, eta, z_part, "WBm2", WBm2);
  } else {
    printf("Weird functional encountered:  %s\n", argv[4]);
    return 1;
  }

  return 0;
}
