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


void run_spherical_solute(double diam, double eta, const char *name, Functional fhs) {
  diameter = diam;
  printf("Diameter is %g hard spheres, filling fraction %g\n", diameter, eta);

  const double myvolume = M_PI*(diameter+1)*(diameter+1)*(diameter+1)/6;
  const double meandensity = eta/(4*M_PI/3);

  Functional f = OfEffectivePotential(HardSpheresNoTensor(1.0) + IdealGas());
  double mu = find_chemical_potential(f, 1, meandensity);
  f = OfEffectivePotential(HardSpheresNoTensor(1.0) + IdealGas()
                           + ChemicalPotential(mu));

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
    double mean_contact_density = ContactDensitySimplest(1.0).integral(1, density)/myvolume;
    
    printf("%g\t%.15g\t%.15g\n", diameter, energy, mean_contact_density);
  }
  
  Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }
  char *plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/contact/figs/inner-sphere%s-%04.1f-%04.2f.dat", name, diameter, eta);
  printf("Saving as %s\n", plotname);
  Grid energy_density(gd, f(1, gd, potential));
  Grid contact_density(gd, ContactDensitySimplest(1.0)(1, gd, density));
  Grid contact_density_sphere(gd, ContactDensitySphere(1.0)(1, gd, density));
  if (strlen(name) == 4) contact_density_sphere = ContactDensitySphereWBm2(1.0)(1, gd, density);
  Grid gross_density(gd, GrossContactDensity(1.0)(1, gd, density));
  Grid n0(gd, ShellConvolve(1)(1, density)/(4*M_PI));
  Grid wu_contact_density(gd, FuWuContactDensity(1.0)(1, gd, density));
  //Grid wu_contact_density_no_zeta(gd, FuWuContactDensityNoZeta(1.0)(1, gd, density));
  radial_plot(plotname, density, energy_density, contact_density, wu_contact_density,
              contact_density_sphere, n0, gross_density);
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
  if (strlen(argv[3]) == 2) {
    run_spherical_solute(diameter, eta, "WB", WB);
  } else if (strlen(argv[3]) == 3) {
    run_spherical_solute(diameter, eta, "WBT", WBT);
  } else if (strlen(argv[3]) == 4) {
    run_spherical_solute(diameter, eta, "WBm2", WBm2);
  } else {
    printf("Weird functional encountered:  %s\n", argv[3]);
    return 1;
  }

  return 0;
}
