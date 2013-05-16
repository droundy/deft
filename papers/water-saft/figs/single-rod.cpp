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
#include "utilities.h"
#include "handymath.h"

const double nmtobohr = 18.8972613; // Converts nm to bohr
const double nm = 18.8972613;
// Here we set up the lattice.
//const double zmax = 5*nm;
//const double ymax = 5*nm;
const double width = 0.0001;
double diameter = 1.0*nm;
double dr = 0.04*nm;

double notinwall(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  if (sqrt(sqr(z)+sqr(y)) < diameter/2) {
      return 0; 
  }
  return 1;
}

void plot_grids_y_direction(const char *fname, const Grid &a) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = gd.Nx/2;
  const int z = 0;
  for (int y=0; y<gd.Ny/2; y++) {
    Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
    double ahere = a(x,y,z);
    fprintf(out, "%g\t%g\t%g\t%g\n", here[0], here[1], here[2], ahere);
  }
  fclose(out);
}

void plot_grids_yz_directions(const char *fname, const Grid &a, const Grid &b, 
			    const Grid &c, const Grid &d) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = gd.Nx/2;
  //const int y = gd.Ny/2;
  for (int y=-gd.Ny/2; y<=gd.Ny/2; y++) {
    for (int z=-gd.Nz/2; z<=gd.Nz/2; z++) {
      Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
      double ahere = a(here);
      double bhere = b(here);
      double chere = c(here);
      double dhere = d(here);
      fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[0], here[1], here[2], 
	      ahere, bhere, chere, dhere);
    }
    fprintf(out,"\n");
 }  
  fclose(out);
}

int main(int argc, char **argv) {
  if (argc == 2) {
    if (sscanf(argv[1], "%lg", &diameter) != 1) {
      printf("Got bad argument: %s\n", argv[1]);
      return 1;
    }
    diameter *= nm;
  } else {
    printf("Usage: %s RADIUS (got argc of %d)\n", argv[0], argc);
    return 1;
  }
  //printf("Diameter is %g bohr = %g nm\n", diameter, diameter/nm);
  const double padding = 2*nm;
  const double ymax = diameter+2*padding;
  const double zmax = diameter+2*padding;

  Functional f = OfEffectivePotential(WaterSaft(new_water_prop.lengthscale,
						new_water_prop.epsilonAB, new_water_prop.kappaAB,
						new_water_prop.epsilon_dispersion,
						new_water_prop.lambda_dispersion,
						new_water_prop.length_scaling, 0));
  double n_1atm = pressure_to_density(f, new_water_prop.kT, atmospheric_pressure,
					      0.001, 0.01);

  double mu_satp = find_chemical_potential(f, new_water_prop.kT, n_1atm);

  f = OfEffectivePotential(WaterSaft(new_water_prop.lengthscale,
				     new_water_prop.epsilonAB, new_water_prop.kappaAB,
				     new_water_prop.epsilon_dispersion,
				     new_water_prop.lambda_dispersion,
				     new_water_prop.length_scaling, mu_satp));
  
  const double EperVolume = f(new_water_prop.kT, -new_water_prop.kT*log(n_1atm));
  const double EperNumber = EperVolume/n_1atm;
  const double EperCell = EperVolume*(zmax*ymax - 0.25*M_PI*diameter*diameter)*width;

  Functional X = Xassociation(new_water_prop.lengthscale, new_water_prop.epsilonAB, 
  			    new_water_prop.kappaAB, new_water_prop.epsilon_dispersion,
  			    new_water_prop.lambda_dispersion,
  			    new_water_prop.length_scaling);
  
  Functional S = OfEffectivePotential(SaftEntropy(new_water_prop.lengthscale, 
						  new_water_prop.epsilonAB, 
						  new_water_prop.kappaAB, 
						  new_water_prop.epsilon_dispersion,
						  new_water_prop.lambda_dispersion,
						  new_water_prop.length_scaling));
    
  Lattice lat(Cartesian(width,0,0), Cartesian(0,ymax,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, 0.1);
    
  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
    
  f = OfEffectivePotential(WaterSaft(new_water_prop.lengthscale,
                                     new_water_prop.epsilonAB, new_water_prop.kappaAB,
                                     new_water_prop.epsilon_dispersion,
                                     new_water_prop.lambda_dispersion,
                                     new_water_prop.length_scaling, mu_satp));
  f = constrain(constraint, f);
  // constraint.epsNativeSlice("papers/water-saft/figs/single-rod-in-water-constraint.eps",
  // 			      Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
  // 			      Cartesian(0,ymax/2,zmax/2));
  //printf("Constraint has become a graph!\n");
   
  potential = new_water_prop.liquid_density*constraint
    + 200*new_water_prop.vapor_density*VectorXd::Ones(gd.NxNyNz);
  //potential = new_water_prop.liquid_density*VectorXd::Ones(gd.NxNyNz);
  potential = -new_water_prop.kT*potential.cwise().log();
  
  const double surface_tension = 5e-5; // crude guess from memory...
  const double surfprecision = 1e-5*M_PI*diameter*width*surface_tension; // five digits accuracy
  const double bulkprecision = 1e-12*fabs(EperCell); // but there's a limit on our precision for small rods
  const double precision = bulkprecision + surfprecision;
  //printf("Precision limit from surface tension is to %g based on %g and %g\n",
  //       precision, surfprecision, bulkprecision);
  Minimizer min = Precision(precision,
                            PreconditionedConjugateGradient(f, gd, new_water_prop.kT,
                                                            &potential,
                                                            QuadraticLineMinimizer));
    
  //printf("\nDiameter of rod = %g bohr (%g nm), dr = %g nm\n", diameter, diameter/nm, dr/nm);
    
  const int numiters = 200;
  for (int i=0;i<numiters && min.improve_energy(true);i++) {
    fflush(stdout);
    //Grid density(gd, EffectivePotentialToDensity()(new_water_prop.kT, gd, potential));
     
    //density.epsNativeSlice("papers/water-saft/figs/single-rod-in-water.eps", 
    //			     Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
    //			     Cartesian(0,ymax/2,zmax/2));
      
    // sleep(3);
    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
    }
  }

  Grid density(gd, EffectivePotentialToDensity()(new_water_prop.kT, gd, potential));
  //printf("The bulk energy per cell should be %g\n", EperCell);
  //printf("The bulk energy based on number should be %g\n", EperNumber*density.integrate());
  //printf("Number of water molecules is %g\n", density.integrate());
  double energy = (min.energy() - EperCell)/width;
  energy = (min.energy() - EperNumber*density.integrate())/width;
  //printf("Energy is %.15g\n", energy);

  char *datname = new char[1024];
  sprintf(datname, "papers/water-saft/figs/single-rod-%04.2fnm-energy.dat", diameter/nm);
  FILE *o = fopen(datname, "w");
  delete[] datname;
  fprintf(o, "%g\t%.15g\n", diameter/nm, energy);
  fclose(o);

  {
    //double peak = peak_memory()/1024.0/1024;
    //double current = current_memory()/1024.0/1024;
    //printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }

  char *plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/water-saft/figs/single-rod-slice-%04.1f.dat", diameter/nm);
  plot_grids_y_direction(plotname, density);
  free(plotname);

  {
    //double peak = peak_memory()/1024.0/1024;
    //double current = current_memory()/1024.0/1024;
    //printf("Peak memory use is %g M (current is %g M)\n", peak, current);
  }
}
