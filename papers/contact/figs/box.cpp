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
#include "Functionals.h"
#include "equation-of-state.h"
#include "LineMinimizer.h"
#include "ContactDensity.h"
#include "utilities.h"
#include "handymath.h"

const double nm = 18.8972613;
// Here we set up the lattice.
double xmax = 3, ymax = 3, zmax = 3;
double N = 4;
bool using_default_box = true;

//Functional HS = HardSpheresNoTensor(1.0);
Functional HS = HardSpheresFast(1.0);

const int numiters = 25;

double N_from_mu(Minimizer *min, Grid *potential, const Grid &constraint, double mu) {
  Functional f = constrain(constraint, OfEffectivePotential(HS + IdealGas()
                                                            + ChemicalPotential(mu)));
  double Nnow = 0;
  min->minimize(f, potential->description());
  for (int i=0;i<numiters && min->improve_energy(false);i++) {
    Grid density(potential->description(), EffectivePotentialToDensity()(1, potential->description(), *potential));
    Nnow = density.sum()*potential->description().dvolume;
    printf("Nnow is %g vs %g\n", Nnow, N);
    fflush(stdout);
    
    density.epsNativeSlice("papers/contact/figs/box.eps", 
                           Cartesian(0,ymax+2,0), Cartesian(0,0,zmax+2), 
                           Cartesian(0,-ymax/2-1,-zmax/2-1));
    density.epsNativeSlice("papers/contact/figs/box-diagonal.eps", 
                           Cartesian(xmax+2,0,zmax+2),  Cartesian(0,ymax+2,0),
                           Cartesian(-xmax/2-1,-ymax/2-1,-zmax/2-1));
    //sleep(3);
  }
  return Nnow;
}

double notinwall(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  const double x = r.x();
  if (fabs(z) < zmax/2 && fabs(y) < ymax/2 && fabs(x) < xmax/2) {
      return 1;
  }
  return 0;
}

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  double peak = peak_memory()/1024.0/1024;
  printf("\t\t%s took %g seconds (peak memory: %g M)\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  last_time = t;
}

void x_plot(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int z = 0;
  const int y = 0;
  for (int x=-gd.Nx/2; x<gd.Nx/2; x++) {
    Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
    double ahere = a(here);
    double bhere = b(here);
    double chere = c(here);
    double dhere = d(here);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\n", here[0], ahere, bhere, chere, dhere);
  }
  fclose(out);
}

void plot_grids_100_center(const char *fname, const Grid &a, const Grid &b, const Grid &c) {
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

void plot_grids_110(const char *fname, const Grid &a, const Grid &b, const Grid &c) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  for (int y=-gd.Ny/2; y<=gd.Ny/2; y++) {
    for (int z=-gd.Nz/2; z<=gd.Nz/2; z++) {
      const int x = y;
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

void plot_grids_100_side(const char *fname, const Grid &a, const Grid &b, 
			    const Grid &c) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  for (int y=-gd.Ny/2; y<=gd.Ny/2; y++) {
    for (int z=-gd.Nz/2; z<=gd.Nz/2; z++) {
      Cartesian here = gd.fineLat.toCartesian(Relative(0,y,z));
      here[0] = xmax/2 - gd.dx;
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

int main(int argc, char *argv[]) {
  if (argc == 5) {
    if (sscanf(argv[1], "%lg", &xmax) != 1) {
      printf("Got bad x argument: %s\n", argv[1]);
      return 1;
    }
    if (sscanf(argv[2], "%lg", &ymax) != 1) {
      printf("Got bad y argument: %s\n", argv[2]);
      return 1;
    }
    if (sscanf(argv[3], "%lg", &zmax) != 1) {
      printf("Got bad z argument: %s\n", argv[3]);
      return 1;
    }
    if (sscanf(argv[4], "%lg", &N) != 1) {
      printf("Got bad N argument: %s\n", argv[4]);
      return 1;
    }
    using_default_box = false;
    printf("Box is %g x %g x %g hard sphere diameters, and it holds %g of them\n", xmax, ymax, zmax, N);
  }

  char *datname = (char *)malloc(1024);
  sprintf(datname, "papers/contact/figs/box-%02.0f,%02.0f,%02.0f-%02.0f-energy.dat", xmax, ymax, zmax, N);
  
  FILE *o = fopen(datname, "w");

  const double myvolume = xmax*ymax*zmax;
  const double meandensity = N/myvolume;

  Functional f = OfEffectivePotential(HS + IdealGas());
  double mu = find_chemical_potential(f, 1, meandensity);
  f = OfEffectivePotential(HS + IdealGas()
                           + ChemicalPotential(mu));

  Lattice lat(Cartesian(xmax+3,0,0), Cartesian(0,ymax+3,0), Cartesian(0,0,zmax+3));
  GridDescription gd(lat, 0.05);
    
  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
  took("Setting the constraint");

  printf("xmax = %g\nymax = %g\nzmax = %g\nmeandensity=%g\n", xmax, ymax, zmax, meandensity);
  f = constrain(constraint, f);
  constraint.epsNativeSlice("papers/contact/figs/box-constraint.eps",
   			      Cartesian(0,ymax+4,0), Cartesian(0,0,zmax+4), 
   			      Cartesian(0,-ymax/2-2,-zmax/2-2));
  printf("Constraint has become a graph!\n");
  
  potential = meandensity*constraint + 1e-4*meandensity*VectorXd::Ones(gd.NxNyNz);
  potential = -potential.cwise().log();
    
  Minimizer min = Precision(1e-6, 
                            PreconditionedConjugateGradient(f, gd, 1, 
                                                            &potential,
                                                            QuadraticLineMinimizer));
    
  double mumax = mu, mumin = mu, dmu = 4.0/N;
  double Nnow = N_from_mu(&min, &potential, constraint, mu);
  if (Nnow > N) {
    while (Nnow > N) {
      mumin = mumax;
      mumax += dmu;
      dmu *= 2;

      Nnow = N_from_mu(&min, &potential, constraint, mumax);
      // Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
      // density = EffectivePotentialToDensity()(1, gd, potential);
      // density.epsNativeSlice("papers/contact/figs/box.eps", 
      //                        Cartesian(0,ymax+2,0), Cartesian(0,0,zmax+2), 
      //                        Cartesian(0,-ymax/2-1,-zmax/2-1));
      // density.epsNativeSlice("papers/contact/figs/box-diagonal.eps", 
      //                        Cartesian(xmax+2,0,zmax+2),  Cartesian(0,ymax+2,0),
      //                        Cartesian(-xmax/2-1,-ymax/2-1,-zmax/2-1));
      printf("mumax %g gives N %g\n", mumax, Nnow);
      took("Finding N from mu");
    }
    printf("mu is between %g and %g\n", mumin, mumax);
  } else {
    while (Nnow < N) {
      mumax = mumin;
      if (mumin > dmu) {
        mumin -= dmu;
        dmu *= 2;
      } else if (mumin > 0) {
        mumin = -mumin;
      } else {
        mumin *= 2;
      }

      Nnow = N_from_mu(&min, &potential, constraint, mumin);
      // density = EffectivePotentialToDensity()(1, gd, potential);
      // density.epsNativeSlice("papers/contact/figs/box.eps", 
      //                        Cartesian(0,ymax+2,0), Cartesian(0,0,zmax+2), 
      //                        Cartesian(0,-ymax/2-1,-zmax/2-1));
      // density.epsNativeSlice("papers/contact/figs/box-diagonal.eps", 
      //                        Cartesian(xmax+2,0,zmax+2),  Cartesian(0,ymax+2,0),
      //                        Cartesian(-xmax/2-1,-ymax/2-1,-zmax/2-1));
      printf("mumin %g gives N %g\n", mumin, Nnow);
      took("Finding N from mu");
    }
    printf("mu is between %g and %g\n", mumin, mumax);
  }

  while (fabs(N/Nnow-1) > 1e-3) {
    mu = 0.5*(mumin + mumax);
    Nnow = N_from_mu(&min, &potential, constraint, mu);
    // density = EffectivePotentialToDensity()(1, gd, potential);
    // density.epsNativeSlice("papers/contact/figs/box.eps", 
    //                        Cartesian(0,ymax+2,0), Cartesian(0,0,zmax+2), 
    //                        Cartesian(0,-ymax/2-1,-zmax/2-1));
    // density.epsNativeSlice("papers/contact/figs/box-diagonal.eps", 
    //                        Cartesian(xmax+2,0,zmax+2),  Cartesian(0,ymax+2,0),
    //                        Cartesian(-xmax/2-1,-ymax/2-1,-zmax/2-1));
    printf("Nnow is %g vs %g with mu %g\n", Nnow, N, mu);
    took("Finding N from mu");
    if (Nnow > N) {
      mumin = mu;
    } else {
      mumax = mu;
    }
  }
  printf("N final is %g (vs %g) with mu = %g\n", Nnow, N, mu);

  double energy = min.energy();
  printf("Energy is %.15g\n", energy);

  Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
  double mean_contact_density = ContactDensitySimplest(1.0).integral(1, density)/myvolume;
  
  fprintf(o, "%g\t%g\t%g\t%.15g\t%.15g\n", xmax, ymax, zmax, energy, mean_contact_density);
  
  Grid energy_density(gd, f(1, gd, potential));
  Grid contact_density(gd, ContactDensitySimplest(1.0)(1, gd, density));
  Grid n0(gd, ShellConvolve(1)(1, density));
  Grid wu_contact_density(gd, FuWuContactDensity(1.0)(1, gd, density));
  char *plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/contact/figs/box-100c--%02.0f,%02.0f,%02.0f-%02.0f.dat", xmax, ymax, zmax, N);
  plot_grids_100_center(plotname, density, energy_density, contact_density);
  sprintf(plotname, "papers/contact/figs/box-100s--%02.0f,%02.0f,%02.0f-%02.0f.dat", xmax, ymax, zmax, N);
  plot_grids_100_side(plotname, density, energy_density, contact_density);
  sprintf(plotname, "papers/contact/figs/box-110c--%02.0f,%02.0f,%02.0f-%02.0f.dat", xmax, ymax, zmax, N);
  plot_grids_110(plotname, density, energy_density, contact_density);
  sprintf(plotname, "papers/contact/figs/box-x-%02.0f,%02.0f,%02.0f-%02.0f.dat", xmax, ymax, zmax, N);
  x_plot(plotname, density, energy_density, contact_density, wu_contact_density);
  free(plotname);
  density.epsNativeSlice("papers/contact/figs/box.eps", 
                         Cartesian(0,ymax+2,0), Cartesian(0,0,zmax+2), 
                         Cartesian(0,-ymax/2-1,-zmax/2-1));
  density.epsNativeSlice("papers/contact/figs/box-diagonal.eps", 
                         Cartesian(xmax+2,0,zmax+2),  Cartesian(0,ymax+2,0),
                         Cartesian(-xmax/2-1,-ymax/2-1,-zmax/2-1));
  
  took("Plotting stuff");
  
  fclose(o);
}
