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

Functional SoftFluid(double radius, double V0, double mu);
Functional HardFluid(double radius, double mu);

// Here we set up the lattice.
static double width = 15;
const double dx = 0.001;
const double dw = 0.001;
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
    fprintf(out, "%g\t%g\n", here[2], ahere);
  }
  fclose(out);
}

double run_walls(double eta, const char *name, Functional fhs, double teff) {
  //printf("Filling fraction is %g with functional %s at temperature %g\n", eta, name, teff);
  //fflush(stdout);
  double kT = teff;
  if (kT == 0) kT = 1;

  Functional f = OfEffectivePotential(fhs);

  const double zmax = width + 2*spacing;
  Lattice lat(Cartesian(dw,0,0), Cartesian(0,dw,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, dx);

  Grid constraint(gd);
  constraint.Set(notinwall);
  f = constrain(constraint, f);

  // We reuse the potential, which should give us a better starting
  // guess on each calculation.
  static Grid *potential = 0;
  static double old_temperature = 0;
  if (strcmp(name, "hard") == 0) {
    // start over for each potential
    delete potential;
    potential = 0;
  }
  if (!potential) {
    potential = new Grid(gd);
    *potential = (eta*constraint + 1e-4*eta*VectorXd::Ones(gd.NxNyNz))/(4*M_PI/3);
    *potential = -kT*potential->cwise().log();
  } else {
    // Adjust the potential so the initial guess for density is the
    // same as we just found in our last simulation.

    *potential *= kT/old_temperature;
  }

  // FIXME below I use the HS energy because of issues with the actual
  // functional.
  const double approx_energy = fhs(kT, eta/(4*M_PI/3))*dw*dw*width;
  const double precision = fabs(approx_energy*1e-5);
  printf("\tMinimizing to %g absolute precision from %g from %g...\n", precision, approx_energy, kT);
  fflush(stdout);

  Minimizer min = Precision(precision,
                            PreconditionedConjugateGradient(f, gd, kT,
                                                            potential,
                                                            QuadraticLineMinimizer));
  took("Setting up the variables");
  if (strcmp(name, "hard") != 0) {
    printf("For now, SoftFluid doesn't work properly, so we're skipping the\n");
    printf("minimization at temperature %g.\n", teff);
  } else {
    for (int i=0;min.improve_energy(false) && i<100;i++) {
    }
  }
  took("Doing the minimization");
  min.print_info();

  Grid density(gd, EffectivePotentialToDensity()(kT, gd, *potential));
  //printf("# per area is %g at filling fraction %g\n", density.sum()*gd.dvolume/dw/dw, eta);

  char *plotname = (char *)malloc(1024);

  sprintf(plotname, "papers/fuzzy-fmt/figs/walls%s-%06.4f-%04.2f.dat", name, teff, eta);
  z_plot(plotname, Grid(gd, 4*M_PI*density/3));
  free(plotname);

  {
    GridDescription gdj = density.description(); 
    double sep =  gdj.dz*gdj.Lat.a3().norm();
    int div = gdj.Nz;
    int mid = int (div/2.0);
    double Ntot_per_A = 0;
    double mydist = 0;
   
    for (int j=0; j<mid; j++){
      Ntot_per_A += density(0,0,j)*sep;
      mydist += sep;
    }

    double Extra_per_A = Ntot_per_A - eta/(4.0/3.0*M_PI)*width/2;

    FILE *fout = fopen("papers/fuzzy-fmt/figs/wallsfillingfracInfo.txt", "a");
    fprintf(fout, "walls%s-%04.2f.dat  -  If you want to match the bulk filling fraction of figs/walls%s-%04.2f.dat, than the number of extra spheres per area to add is %04.10f.  So you'll want to multiply %04.2f by your cavity volume and divide by (4/3)pi.  Then add %04.10f times the Area of your cavity to this number\n",
	    name, eta, name, eta, Extra_per_A, eta, Extra_per_A);

    int wallslen = 20;
    double Extra_spheres =  (eta*wallslen*wallslen*wallslen/(4*M_PI/3) + Extra_per_A*wallslen*wallslen);  
    fprintf (fout, "For filling fraction %04.02f and walls of length %d you'll want to use %.0f spheres.\n\n", eta, wallslen, Extra_spheres);

    fclose(fout); 
  }

  {
    //double peak = peak_memory()/1024.0/1024;
    //double current = current_memory()/1024.0/1024;
    //printf("Peak memory use is %g M (current is %g M)\n", peak, current);

  }

  old_temperature = kT;

  took("Plotting stuff");
  printf("density %g gives ff %g for eta = %g and T = %g\n", density(0,0,gd.Nz/2),
         density(0,0,gd.Nz/2)*4*M_PI/3, eta, teff);
  return density(0, 0, gd.Nz/2)*4*M_PI/3; // return bulk filling fraction
}

int main(int, char **) {
  FILE *fout = fopen("papers/fuzzy-fmt/figs/wallsfillingfracInfo.txt", "w");
  fclose(fout);
  const double etas[] = { 0.1, 0.2, 0.4 };
  const double temps[] = { 0.0, 0.01, 0.02, 0.03 };
  for (double eta = 0.5; eta > 0; eta-=0.1) {
    for (unsigned int i = 0; i<sizeof(temps)/sizeof(temps[0]); i++) {
      const double temp = temps[i];
      Functional f = HardFluid(1,0);
      if (temp > 0) f = SoftFluid(1, 1, 0);
      const double mu = find_chemical_potential(OfEffectivePotential(f), (temp)?temp:1, eta/(4*M_PI/3));
      printf("mu is %g for eta = %g at temperature %g\n", mu, eta, temp);
      if (temp > 0) f = SoftFluid(1, 1, mu);
      else f = HardFluid(1, mu);

      const char *name = "hard";
      if (temp > 0) name = "soft";

      run_walls(eta, name, f, temp);
    }
  }
  // Just create this file so make knows we have run.
  if (!fopen("papers/fuzzy-fmt/figs/walls.dat", "w")) {
    printf("Error creating walls.dat!\n");
    return 1;
  }
  return 0;
}
