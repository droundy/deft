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
static double width = 30;
const double dw = 0.001;
const double spacing = 3; // space on each side

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
  last_time = t;
}

Functional WB = HardSpheresNoTensor(1.0);
Functional WBm2 = HardSpheresWBm2(1.0);
Functional WBT = HardSpheresWBFast(1.0);

const int numiters = 25;

void z_plot(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d,
            const Grid &e, const Grid &f, const Grid &g) {
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
    double bhere = b(x,y,z);
    double chere = c(x,y,z);
    double dhere = d(x,y,z);
    double ehere = e(x,y,z);
    double fhere = f(x,y,z);
    double ghere = g(x,y,z);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[2],
            ahere, bhere, chere, dhere, ehere, fhere, ghere);
  }
  fclose(out);
}

void plot_grids_yz_directions(const char *fname, const Grid &a, const Grid &b, 
			    const Grid &c) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = 0;
  //const int y = gd.Ny/2;
  for (int y=-gd.Ny/2; y<=gd.Ny/2; y++) {
    for (int z=-gd.Nz/2; z<=gd.Nz/2; z++) {
      Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
      double ahere = a(here);
      double bhere = b(here);
      double chere = c(here);
      fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\n", here[0], here[1], here[2], 
	      ahere, bhere, chere);
    }
    fprintf(out,"\n");
 }  
  fclose(out);
}


void run_walls(double eta, const char *name, Functional fhs) {
  printf("Filling fraction is %g\n", eta);

  Functional f = OfEffectivePotential(fhs + IdealGas());
  double mu = find_chemical_potential(f, 1, eta/(4*M_PI/3));
  f = OfEffectivePotential(fhs + IdealGas()
                           + ChemicalPotential(mu));

  const double zmax = width + 2*spacing;
  Lattice lat(Cartesian(dw,0,0), Cartesian(0,dw,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, 0.01);
    
  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
  f = constrain(constraint, f);
  
  potential = (eta*constraint + 1e-4*eta*VectorXd::Ones(gd.NxNyNz))/(4*M_PI/3);
  potential = -potential.cwise().log();
  //f.run_finite_difference_test("foobar", 1, potential);

  const double approx_energy = (fhs + IdealGas() + ChemicalPotential(mu))(1, eta/(4*M_PI/3))*dw*dw*width;
  const double precision = fabs(approx_energy*1e-4);
  printf("Minimizing to %g absolute precision...\n", precision);

  Minimizer min = Precision(precision,
                            PreconditionedConjugateGradient(f, gd, 1, 
                                                            &potential,
                                                            QuadraticLineMinimizer));
  for (int i=0;min.improve_energy(false) && i<100;i++) {
  }
  took("Doing the minimization");
    
  double energy = min.energy();
  printf("Energy is %.15g\n", energy);
  
  Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
  double mean_contact_density = ContactDensitySimplest(1.0).integral(1, density)/(dw*dw*width);
    
  printf("%g\t%.15g\t%.15g\n", width, energy, mean_contact_density);
  
  char *plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/contact/figs/walls%s-%04.1f-%04.2f.dat", name, width, eta);
  Grid energy_density(gd, f(1, gd, potential));
  Grid contact_density_S(gd, ContactDensity_S(1.0)(1, gd, density));
  Grid contact_density_sphere(gd, ContactDensitySphere(1.0)(1, gd, density));
  if (strlen(name) == 4) contact_density_sphere = ContactDensitySphereWBm2(1.0)(1, gd, density);
  Grid gross_density(gd, GrossContactDensity(1.0)(1, gd, density));
  Grid n0(gd, ShellConvolve(1)(1, density)/(4*M_PI));
  Grid wu_contact_density(gd, YuWuContact(1.0)(1, gd, density));
  //Grid wu_contact_density_no_zeta(gd, FuWuContactDensityNoZeta(1.0)(1, gd, density));
  // plot_grids_yz_directions(plotname, density, energy_density, contact_density);
  sprintf(plotname, "papers/contact/figs/walls%s-%04.2f.dat", name, eta);
  z_plot(plotname, density, energy_density, contact_density_S, wu_contact_density, contact_density_sphere,
         n0, gross_density);
  free(plotname);
  // density.epsNativeSlice("papers/contact/figs/walls.eps", 
  //                        Cartesian(0,xmax,0), Cartesian(0,0,xmax), 
  //                        Cartesian(0,xmax/2,xmax/2));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }
  
  took("Plotting stuff");
}

int main(int, char **) {
  for (double eta = 0.1; eta < 0.6; eta+=0.1) {
    run_walls(eta, "WB", WB);
    run_walls(eta, "WBT", WBT);
    run_walls(eta, "WBm2", WBm2);
  }
  // Just create this file so make knows we have run.
  if (!fopen("papers/contact/figs/walls.dat", "w")) {
    printf("Error creating walls.dat!\n");
    return 1;
  }
  return 0;
}
