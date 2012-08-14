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
#include "OptimizedFunctionals.h"
#include "equation-of-state.h"
#include "LineMinimizer.h"
#include "ContactDensity.h"
#include "utilities.h"
#include "handymath.h"

// Here we set up the lattice.
static double diameter;
const double padding = 16;

double notinwall(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  const double x = r.x();
  if (sqrt(sqr(z)+sqr(y)+sqr(x)) > diameter/2) {
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
        if ((here-pos-delta).dot(delta) > 0 && (here-pos-delta).norm() < closest) {
          //printf("here we are... %d %d %d\n", x, y, z);
          bestyet = g(Cartesian(pos + here));
          closest = (here-pos-delta).norm();
        }
        //printf("Got here now..;\n");
      }
    }
  }
  return bestyet;
}
// The following computes the average value for the grid x at radius
// from pos.
double radial_average(const Grid &g, Cartesian pos, double radius) {
  double total = 0;
  for (int i=0;i<20;i++) {
    total += closest_value(g, pos, Vector3d(radius * icosohedron[i]/icosohedron[i].norm()));
  }
  return total/20;
}

void radial_plot2(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d, const Grid &e,
                 const Grid &f, const Grid &g) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const double dr = diameter / floor(diameter/gd.fineLat.a2().norm());
  //fprintf(out, "dr = %g\n", dr);
  //fprintf(out, "gd.dx = %g\n", gd.dx);
  //fprintf(out, "Nx*dr = %g\n", gd.Nx*dr);
  //fprintf(out, "Ny*dr = %g\n", gd.Ny*dr);
  for (double r = 0; r < gd.Nx*dr/2; r += dr) {
    //printf("WOrking on %g\n", r);
    const double meana = radial_average(a, Cartesian(0,0,0), r);
    const double meanb = radial_average(b, Cartesian(0,0,0), r);
    const double meanc = radial_average(c, Cartesian(0,0,0), r);
    const double meand = radial_average(d, Cartesian(0,0,0), r);
    const double meane = radial_average(e, Cartesian(0,0,0), r);
    const double meanf = radial_average(f, Cartesian(0,0,0), r);
    const double meang = radial_average(g, Cartesian(0,0,0), r);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", r,
            meana,
            meanb,
            meanc,
            meand,
            meane,
            meanf,
            meang);
  }
  fclose(out);
}

void radial_plot(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d, const Grid &e,
                 const Grid &f, const Grid &g) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = 0;
  const int z = 0;
  for (int y=0; y<gd.Ny/2; y++) {
    Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
    printf("r = %g vs %g\n", here.norm(), y*gd.dx);
    double ahere = a(x,y,z);
    double bhere = b(x,y,z);
    double chere = c(x,y,z);
    double dhere = d(x,y,z);
    double ehere = e(x,y,z);
    double fhere = f(x,y,z);
    double ghere = g(x,y,z);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[1],
            ahere, bhere, chere, dhere, ehere, fhere, ghere);
  }
  fclose(out);
}


void run_spherical_solute(double diam, double eta, const char *name, Functional fhs) {
  diameter = diam;
  printf("Diameter is %g hard spheres, filling fraction %g\n", diameter, eta);

  const double meandensity = eta/(4*M_PI/3);

  Functional f = OfEffectivePotential(fhs + IdealGas());
  double mu = find_chemical_potential(f, 1, meandensity);
  f = OfEffectivePotential(fhs + IdealGas() + ChemicalPotential(mu));

  const double xmax = diameter + padding;
  Lattice lat(Cartesian(xmax,0,0), Cartesian(0,xmax,0), Cartesian(0,0,xmax));
  GridDescription gd(lat, 0.1);

  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
  took("Setting the constraint");

  f = constrain(constraint, f);

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
  printf("N = %g\n", density.sum()*gd.dvolume);
  char *plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/contact/figs/inner-sphere%s-%04.1f-%04.2f.dat", name, diameter, eta);
  printf("Saving as %s\n", plotname);
  Grid correlation_S(gd, Correlation_S2(1.0)(1, gd, density));
  Grid correlation_A(gd, Correlation_A2(1.0)(1, gd, density));
  if (strlen(name) == 4) {
    printf("Computing correlation for mark II version...\n");
    correlation_S = Correlation_S_WBm2(1.0)(1, gd, density);
    correlation_A = Correlation_A_WBm2(1.0)(1, gd, density);
  }
  Grid gross_correlation(gd, GrossCorrelation(1.0)(1, gd, density));
  Grid n0(gd, ShellConvolve(1)(1, density)/(4*M_PI));
  Grid nA(gd, ShellConvolve(2)(1, density)/(4*M_PI*4));
  Grid yuwu_correlation(gd, YuWuCorrelation_S(1.0)(1, gd, density));
  sprintf(plotname, "papers/contact/figs/inner-sphere%s-%04.1f-%04.2f.dat", name, diameter, eta);
  radial_plot(plotname, density, n0, correlation_S, yuwu_correlation,
              nA, correlation_A, gross_correlation);
  sprintf(plotname, "papers/contact/figs/inner-sphere%s-%04.1f-%04.2f-mean.dat", name, diameter, eta);
  radial_plot2(plotname, density, n0, correlation_S, yuwu_correlation,
               nA, correlation_A, gross_correlation);
  free(plotname);

  {
    const GridDescription gdp = density.description();
    double inner_rad = diameter/2.0;

    double Ntot = density.sum()*gdp.dvolume;
    double Ndisplaced = eta*gdp.Lat.volume()/(4*M_PI/3) - Ntot;

    double mc_side_len = 25;
    double N = eta*mc_side_len*mc_side_len*mc_side_len/(4.0/3.0*M_PI) - Ndisplaced;

    FILE *fout = fopen("papers/contact/figs/innerfillingfracInfo.txt", "a");
    if (fout==0){
      printf("Not able to open papers/contact/figs/innerfillingfracInfo.txt properly\n");
      exit (1);
    }

    fprintf (fout, "For filling fraction %04.02f, inner-sphere size %04.02f and walls of length %04.02f you'll want to use %.0f spheres.\n\n", eta, inner_rad, mc_side_len, N);

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
  if (argc != 4) {
    printf("got argc %d\n", argc);
    printf("usage: %s diameter filling_fraction (WB|WBT|WBm2)\n", argv[0]);
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
  FILE *fout = fopen("papers/contact/figs/innerfillingfracInfo.txt", "w");
  if (!fout) {
    printf("Unable to create innerfillingfracInfo.txt!\n");
    exit(1);
  }
  fclose(fout);
  if (strlen(argv[3]) == 2) {
    printf("Using Mark I version of White Bear functional...\n");
    run_spherical_solute(diameter, eta, "WB", WB);
  } else if (strlen(argv[3]) == 3) {
    printf("Using Mark I tensorized version of White Bear functional...\n");
    run_spherical_solute(diameter, eta, "WBT", WBT);
  } else if (strlen(argv[3]) == 4) {
    printf("Using Mark II version of White Bear functional...\n");
    run_spherical_solute(diameter, eta, "WBm2", WBm2);
  } else {
    printf("Weird functional encountered:  %s\n", argv[3]);
    return 1;
  }

  return 0;
}
