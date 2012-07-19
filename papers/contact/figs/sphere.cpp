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

const double nm = 18.8972613;
// Here we set up the lattice.
static double cavity_radius;
const double padding = 4;

double notinwall(Cartesian r) {
  if (r.norm() < cavity_radius) {
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

Functional WB = HardSpheresNoTensor2(1.0);
Functional WBT = HardSpheresWBFast(1.0);
Functional WBm2 = HardSpheresWBm2(1.0);

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


void run_spherical_solute(double rad, double eta, const char *name, Functional fhs) {
  cavity_radius = rad;
  printf("Cavity radius is %g hard spheres, filling fraction %g\n", cavity_radius, eta);

  const double meandensity = eta/(4*M_PI/3);

  Functional f = OfEffectivePotential(fhs + IdealGas());
  double mu = find_chemical_potential(f, 1, meandensity);
  f = OfEffectivePotential(fhs + IdealGas() + ChemicalPotential(mu));

  const double xmax = 2*cavity_radius + padding;
  Lattice lat(Cartesian(xmax,0,0), Cartesian(0,xmax,0), Cartesian(0,0,xmax));
  GridDescription gd(lat, 0.1);

  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
  took("Setting the constraint");
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
    fflush(stdout);
  }

  f = constrain(constraint, f);

  potential = meandensity*constraint + 1e-4*meandensity*VectorXd::Ones(gd.NxNyNz);
  potential = -potential.cwise().log();
  //f.run_finite_difference_test("foobar", 1, potential);

  {
    const double approx_energy = (fhs + IdealGas() + ChemicalPotential(mu))(1, eta/(4*M_PI/3))*4*M_PI*rad*rad*rad/3;
    const double precision = fabs(approx_energy*1e-10);
    printf("Minimizing to %g absolute precision...\n", precision);
    Minimizer min = Precision(precision,
                              PreconditionedConjugateGradient(f, gd, 1,
                                                              &potential,
                                                              QuadraticLineMinimizer));
    for (int i=0;min.improve_energy(true) && i<200;i++) {
      took("One iteration");
      {
        double peak = peak_memory()/1024.0/1024;
        double current = current_memory()/1024.0/1024;
        printf("Peak memory use is %g M (current is %g M)\n", peak, current);
        fflush(stdout);
      }
    }
    printf("All done minimizing...");
    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
    }
    fflush(stdout);
    constraint.resize(0); // free memory used by constraint
    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Freed constraint and memory use is %g M (current is %g M)\n", peak, current);
    }

    double energy = min.energy();
    printf("Energy is %.15g\n", energy);
    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
    }
    fflush(stdout);

    printf("%g\t%.15g\n", cavity_radius, energy);
    fflush(stdout);
  }
  Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
  potential.resize(0); // free memory used by potential.
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }
  printf("N = %g\n", density.sum()*gd.dvolume);
  char *plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/contact/figs/outer-sphere%s-%02.0f-%04.1f.dat", name, cavity_radius, eta);
  printf("Saving as %s\n", plotname);
  Grid correlation_S(gd, Correlation_S2(1.0)(1, gd, density));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory after correlation_S is %g M (current is %g M)\n", peak, current);
  }
  Grid correlation_A(gd, Correlation_A2(1.0)(1, gd, density));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory after correlation_A is %g M (current is %g M)\n", peak, current);
  }
  if (strlen(name) == 4) {
    correlation_S = Correlation_S_WBm2(1.0)(1, gd, density);
    correlation_A = Correlation_A_WBm2(1.0)(1, gd, density);
  }
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory after WBM2 is %g M (current is %g M)\n", peak, current);
  }
  Grid gross_correlation(gd, GrossCorrelation(1.0)(1, gd, density));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory after gross_correlation is %g M (current is %g M)\n", peak, current);
  }
  Grid yuwu_correlation(gd, YuWuCorrelation_S(1.0)(1, gd, density));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }
  Grid n0(gd, ShellConvolve(1)(1, density)/(4*M_PI));
  Grid nA(gd, ShellConvolve(2)(1, density)/(4*M_PI*4));
  sprintf(plotname, "papers/contact/figs/outer-sphere%s-%02.0f-%04.1f.dat", name, cavity_radius, eta);
  radial_plot(plotname, density, n0, nA, correlation_S, correlation_A,
              yuwu_correlation, gross_correlation);
  free(plotname);

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
    printf("usage: %s radius filling_fraction (WB|WBT|WBm2)\n", argv[0]);
    return 1;
  }
  if (sscanf(argv[1], "%lg", &cavity_radius) != 1) {
    printf("Got bad cavity_radius argument: %s\n", argv[1]);
    return 1;
  }
  if (sscanf(argv[2], "%lg", &eta) != 1) {
    printf("Got bad eta argument: %s\n", argv[2]);
    return 1;
  }
  if (strlen(argv[3]) == 2) {
    run_spherical_solute(cavity_radius, eta, "WB", WB);
  } else if (strlen(argv[3]) == 3) {
    run_spherical_solute(cavity_radius, eta, "WBT", WBT);
  } else if (strlen(argv[3]) == 4) {
    run_spherical_solute(cavity_radius, eta, "WBm2", WBm2);
  } else {
    printf("Weird functional encountered:  %s\n", argv[3]);
    return 1;
  }

  return 0;
}
