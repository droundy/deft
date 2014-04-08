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

//const double nm = 18.8972613;
// Here we set up the lattice.
double zmax = 25;
double ymax = 25;
double xmax = 25;
double dx = 0.15;
double V0 = 10;
double radius = 1;
double temperature = 0.01; //default temp 
// double kT = 0.01; //default temp 

double soft_sphere_potential(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  const double x = r.x();
  const double distance = sqrt(x*x + y*y + z*z);
  if (distance >= radius*2) return 0;
  return V0*(1-distance/(2*radius))/temperature;
  // return V0*(1-distance/(2*radius))/kT;
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

double run_soft_sphere(double eta, Functional fss, double teff) {
  //printf("Filling fraction is %g with functional %s at temperature %g\n", eta, teff);
  //fflush(stdout);
  temperature = teff;
  //if (kT == 0) kT = ;1

  Functional f = OfEffectivePotential(fss);

  Lattice lat(Cartesian(xmax,0,0), Cartesian(0,ymax,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, dx);

  Grid softspherepotential(gd);
  softspherepotential.Set(soft_sphere_potential);
  
  static Grid *potential = 0;
  potential = new Grid(gd);
  *potential = softspherepotential - temperature*VectorXd::Ones(gd.NxNyNz); // Bad starting guess 
  const double approx_energy = fss(temperature, eta/(4*M_PI/3))*xmax*ymax*zmax;
  const double precision = fabs(approx_energy*1e-8);
  printf("\tMinimizing to %g absolute precision from %g from %g...\n", precision, approx_energy, temperature);
  fflush(stdout);

  Minimizer min = Precision(precision,
                            PreconditionedConjugateGradient(f, gd, temperature,
                                                            potential,
                                                            QuadraticLineMinimizer));
  took("Setting up the variables");
  for (int i=0;min.improve_energy(true) && i<100;i++) {
  }

  took("Doing the minimization");
  min.print_info();

  Grid density(gd, EffectivePotentialToDensity()(temperature, gd, *potential));
  //printf("# per area is %g at filling fraction %g\n", density.sum()*gd.dvolume/dw/dw, eta);

  char *plotname = (char *)malloc(1024);

  sprintf(plotname, "papers/fuzzy-fmt/figs/soft-sphere%06.4f-%04.2f.dat", teff, eta);
  z_plot(plotname, Grid(gd, 4*M_PI*density/3));
  free(plotname);

  // {
  //   GridDescription gdj = density.description(); 
  //   double sep =  gdj.dz*gdj.Lat.a3().norm();
  //   int div = gdj.Nz;
  //   int mid = int (div/2.0);
  //   double Ntot_per_V = 0;
  //   double mydist = 0;
   
  //   for (int j=0; j<mid; j++){
  //     Ntot_per_V += density(0,0,j)*sep;
  //     mydist += sep;
  //   }

  //   double Extra_per_V = Ntot_per_V - eta/(4.0/3.0*M_PI);

  //   FILE *fout = fopen("papers/fuzzy-fmt/figs/wallsfillingfracInfo.txt", "a");
  //   fprintf(fout, "soft-sphere-%04.2f.dat  -  If you want to match the bulk filling fraction of figs/soft-sphere-%04.2f.dat, than the number of extra spheres per area to add is %04.10f.  So you'll want to multiply %04.2f by your cavity volume and divide by (4/3)pi.  Then add %04.10f times the Area of your cavity to this number\n", eta
  //           , eta, Extra_per_V, eta, Extra_per_V);

  //   int wallslen = 20;
  //   double Extra_spheres =  (eta*wallslen*wallslen*wallslen/(4*M_PI/3) + Extra_per_V*wallslen*wallslen);  
  //   fprintf (fout, "For filling fraction %04.02f and walls of length %d you'll want to use %.0f spheres.\n\n", eta, wallslen, Extra_spheres);

  //   fclose(fout); 
  // }

  {
    //double peak = peak_memory()/1024.0/1024;
    //double current = current_memory()/1024.0/1024;
    //printf("Peak memory use is %g M (current is %g M)\n", peak, current);

  }

  took("Plotting stuff");
  printf("density %g gives ff %g for eta = %g and T = %g\n", density(0,0,gd.Nz/2),
         density(0,0,gd.Nz/2)*4*M_PI/3, eta, teff);
  return density(0, 0, gd.Nz/2)*4*M_PI/3; // return bulk filling fraction
}

int main(int argc, char *argv[]) {
  FILE *fout = fopen("papers/fuzzy-fmt/figs/wallsfillingfracInfo.txt", "w");
  fclose(fout);
  const double temps[] = { 0.03, 0.01, 0.001,};
  for (double eta = 0.4; eta >= 0.05; eta-=0.1) {
    for (unsigned int i = 0; i<sizeof(temps)/sizeof(temps[0]); i++) {
      const double temp = temps[i];
      Functional f = SoftFluid(radius, 1, 0);
      const double mu = find_chemical_potential(OfEffectivePotential(f), (temp)?temp:1, eta/(4*M_PI/3));
      printf("mu is %g for eta = %g at temperature %g\n", mu, eta, temp);
      f = SoftFluid(radius, 1, mu);

      run_soft_sphere(eta, f, temp);
    }
  }
  // Just create this file so make knows we have run.
  if (!fopen("papers/fuzzy-fmt/figs/soft-sphere.dat", "w")) {
    printf("Error creating soft-sphere.dat!\n");
    return 1;
  }
  return 0;
}
