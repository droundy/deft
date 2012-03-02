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
static double diameter;

double notinwall(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  const double x = r.x();
  if (sqrt(sqr(z)+sqr(y)+sqr(x)) < diameter/2) {
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

Functional WB = HardSpheresNoTensor(1.0);
Functional WBT = HardSpheresWBFast(1.0);
Functional WBm2 = HardSpheresWBm2(1.0);

const int numiters = 25;

double N_from_mu(Functional fhs, Minimizer *min, Grid *potential,
                 const Grid &constraint, double mu) {
  Functional f = constrain(constraint, OfEffectivePotential(fhs + IdealGas()
                                                            + ChemicalPotential(mu)));
  //double Nnow = 0;
  min->minimize(f, potential->description());
  for (int i=0;i<numiters && min->improve_energy(false);i++) {
    //Grid density(potential->description(), EffectivePotentialToDensity()(1, potential->description(), *potential));
    //Nnow = density.sum()*potential->description().dvolume;
    //printf("Nnow is %g vs %g\n", Nnow, N);
    //fflush(stdout);
    
    //density.epsNativeSlice("papers/contact/figs/sphere.eps", 
    //                       Cartesian(0,xmax,0), Cartesian(0,0,xmax), 
    //                       Cartesian(0,xmax/2,xmax/2));
    //sleep(3);
  }
  Grid density(potential->description(), EffectivePotentialToDensity()(1, potential->description(), *potential));
  return density.sum()*potential->description().dvolume;
}

void radial_plot(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d,
                 const Grid &e, const Grid &f, const Grid &g) {
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


void run_spherical_cavity(double diam, int N, const char *name, Functional fhs) {
  diameter = diam;
  printf("Diameter is %g hard spheres, and it holds %d of them\n", diameter, N);

  char *datname = (char *)malloc(1024);
  sprintf(datname, "papers/contact/figs/sphere%s-%04.1f-%02d-energy.dat", name, diameter, N);
  
  FILE *o = fopen(datname, "w");

  const double myvolume = M_PI*(diameter+1)*(diameter+1)*(diameter+1)/6;
  const double meandensity = N/myvolume;

  Functional f = OfEffectivePotential(HardSpheresNoTensor(1.0) + IdealGas());
  double mu = 2 + find_chemical_potential(f, 1, meandensity);
  f = OfEffectivePotential(HardSpheresNoTensor(1.0) + IdealGas()
                           + ChemicalPotential(mu));

  const double xmax = diameter + 4;
  Lattice lat(Cartesian(xmax,0,0), Cartesian(0,xmax,0), Cartesian(0,0,xmax));
  GridDescription gd(lat, 0.1);
    
  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
  took("Setting the constraint");

  f = constrain(constraint, f);
  //constraint.epsNativeSlice("papers/contact/figs/sphere-constraint.eps",
  // 			      Cartesian(0,xmax,0), Cartesian(0,0,xmax), 
  // 			      Cartesian(0,xmax/2,xmax/2));
  //printf("Constraint has become a graph!\n");
  
  potential = meandensity*constraint + 1e-4*meandensity*VectorXd::Ones(gd.NxNyNz);
  potential = -potential.cwise().log();
  f.run_finite_difference_test("foobar", 1, potential);

  {
    Minimizer min = Precision(1e-6, 
                              PreconditionedConjugateGradient(f, gd, 1, 
                                                              &potential,
                                                              QuadraticLineMinimizer));
    double mumax = mu, mumin = mu, dmu = 4.0/N;
    double Nnow = N_from_mu(fhs, &min, &potential, constraint, mu);
    const double fraccuracy = 1e-3;
    if (fabs(Nnow/N - 1) > fraccuracy) {
      if (Nnow > N) {
        while (Nnow > N) {
          mumin = mumax;
          mumax += dmu;
          dmu *= 2;
          
          Nnow = N_from_mu(fhs, &min, &potential, constraint, mumax);
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
          
          Nnow = N_from_mu(fhs, &min, &potential, constraint, mumin);
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
      
      while (fabs(N/Nnow-1) > fraccuracy) {
        mu = 0.5*(mumin + mumax);
        Nnow = N_from_mu(fhs, &min, &potential, constraint, mu);
        // density = EffectivePotentialToDensity()(1, gd, potential);
        // density.epsNativeSlice("papers/contact/figs/box.eps", 
        //                        Cartesian(0,ymax+2,0), Cartesian(0,0,zmax+2), 
        //                        Cartesian(0,-ymax/2-1,-zmax/2-1));
        // density.epsNativeSlice("papers/contact/figs/box-diagonal.eps", 
        //                        Cartesian(xmax+2,0,zmax+2),  Cartesian(0,ymax+2,0),
        //                        Cartesian(-xmax/2-1,-ymax/2-1,-zmax/2-1));
        printf("Nnow is %g vs %d with mu %g\n", Nnow, N, mu);
        took("Finding N from mu");
        if (Nnow > N) {
          mumin = mu;
        } else {
          mumax = mu;
        }
      }
      printf("N final is %g (vs %d) with mu = %g\n", Nnow, N, mu);
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
    
    fprintf(o, "%g\t%.15g\t%.15g\n", diameter, energy, mean_contact_density);
  }
  
  Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }
  char *plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/contact/figs/sphere%s-%04.1f-%02d.dat", name, diameter, N);
  Grid energy_density(gd, f(1, gd, potential));
  Grid contact_density_S(gd, ContactDensity_S(1.0)(1, gd, density));
  Grid contact_density_sphere(gd, ContactDensitySphere(1.0)(1, gd, density));
  if (strlen(name) == 4) { 
    contact_density_S = ContactDensity_S_WBm2(1.0)(1, gd, density);    
    contact_density_sphere = ContactDensitySphereWBm2(1.0)(1, gd, density);
  }
  Grid gross_density(gd, GrossContactDensity(1.0)(1, gd, density));
  Grid n0(gd, ShellConvolve(1)(1, density)/(4*M_PI));
  Grid wu_contact_density(gd, FuWuContactDensity(1.0)(1, gd, density));
  //Grid wu_contact_density_no_zeta(gd, FuWuContactDensityNoZeta(1.0)(1, gd, density));
  // plot_grids_yz_directions(plotname, density, energy_density, contact_density);
  sprintf(plotname, "papers/contact/figs/sphere%s-radial-%04.1f-%02d.dat", name, diameter, N);
  radial_plot(plotname, density, energy_density, contact_density_S, wu_contact_density,
              contact_density_sphere, n0, gross_density);
  free(plotname);
  // density.epsNativeSlice("papers/contact/figs/sphere.eps", 
  //                        Cartesian(0,xmax,0), Cartesian(0,0,xmax), 
  //                        Cartesian(0,xmax/2,xmax/2));
  {
    double peak = peak_memory()/1024.0/1024;
    double current = current_memory()/1024.0/1024;
    printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }
  
  took("Plotting stuff");
  
  fclose(o);
}

int main(int argc, char *argv[]) {
  int N;
  if (argc != 4) {
    printf("got argc %d\n", argc);
    printf("usage: %s diameter N (WB|WBT|WBm2)\n", argv[0]);
    return 1;
  }
  if (sscanf(argv[1], "%lg", &diameter) != 1) {
    printf("Got bad diameter argument: %s\n", argv[1]);
    return 1;
  }
  if (sscanf(argv[2], "%d", &N) != 1) {
    printf("Got bad N argument: %s\n", argv[2]);
    return 1;
  }
  if (strlen(argv[3]) == 2) {
    run_spherical_cavity(diameter, N, "WB", WB);
  } else if (strlen(argv[3]) == 3) {
    run_spherical_cavity(diameter, N, "WBT", WBT);
  } else if (strlen(argv[3]) == 4) {
    run_spherical_cavity(diameter, N, "WBm2", WBm2);
  } else {
    printf("Weird functional encountered:  %s\n", argv[3]);
    return 1;
  }

  return 0;
}
