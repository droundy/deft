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
#include "utilities.h"
#include "handymath.h"

const double nmtobohr = 18.8972613; // Converts nm to bohr
const double nm = 18.8972613;
// Here we set up the lattice.
const double zmax = 3*nm;
const double ymax = 6*nm;
const double width = 0.0001;
const double cavitysize = 1*nm;
double distance = 1;

double notinwall(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  if (sqrt(sqr(z)+sqr(y+(cavitysize+distance)/2)) < cavitysize/2) {
      return 0; 
  } 
  if (sqrt(sqr(z)+sqr(y-(cavitysize+distance)/2)) < cavitysize/2) {
      return 0;
  }
  return 1;
}

void plot_grids_y_direction(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = gd.Nx/2;
  const int z = 0;
  for (int y=0; y<gd.Ny; y++) {
    Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
    double ahere = a(x,y,z);
    double bhere = b(x,y,z);
    double chere = c(x,y,z);
    double dhere = d(x,y,z);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[0], here[1], here[2], ahere, bhere, chere, dhere);
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
  for (int y=0; y<gd.Ny; y++) {
    for (int z=0; z<gd.Nz; z++) {
      Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
      double ahere = a(x,y,z);
      double bhere = b(x,y,z);
      double chere = c(x,y,z);
      double dhere = d(x,y,z);
      fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[0], here[1], here[2], 
	      ahere, bhere, chere, dhere);
    }
    fprintf(out,"\n");
 }  
  fclose(out);
}

int main(int, char **) {
  FILE *o = fopen("paper/figs/rods-in-water.dat", "w");

  Functional f = OfEffectivePotential(SaftFluid(water_prop.lengthscale,
						water_prop.epsilonAB, water_prop.kappaAB,
						water_prop.epsilon_dispersion,
						water_prop.lambda_dispersion,
						water_prop.length_scaling, 0));
  double n_1atm = pressure_to_density(f, water_prop.kT, atmospheric_pressure,
					      0.001, 0.01);

  double mu_satp = find_chemical_potential(f, water_prop.kT, n_1atm);

  f = OfEffectivePotential(SaftFluid(water_prop.lengthscale,
		water_prop.epsilonAB, water_prop.kappaAB,
		water_prop.epsilon_dispersion,
		water_prop.lambda_dispersion,
				     water_prop.length_scaling, mu_satp));
  
  Functional X = Xassociation(water_prop.lengthscale, water_prop.epsilonAB, 
  			    water_prop.kappaAB, water_prop.epsilon_dispersion,
  			    water_prop.lambda_dispersion,
  			    water_prop.length_scaling);
  
  Functional S = OfEffectivePotential(SaftEntropy(water_prop.lengthscale, 
						  water_prop.epsilonAB, 
						  water_prop.kappaAB, 
						  water_prop.epsilon_dispersion,
						  water_prop.lambda_dispersion,
						  water_prop.length_scaling));
  for (distance=0*nm; distance<2.0*nm; distance += .1*nm) {
    Lattice lat(Cartesian(width,0,0), Cartesian(0,ymax,0), Cartesian(0,0,zmax));
    GridDescription gd(lat, 0.2);
    
    Grid potential(gd);
    Grid constraint(gd);
    constraint.Set(notinwall);
    
    f = OfEffectivePotential(SaftFluid(water_prop.lengthscale,
				       water_prop.epsilonAB, water_prop.kappaAB,
				       water_prop.epsilon_dispersion,
				       water_prop.lambda_dispersion,
				       water_prop.length_scaling, mu_satp));
    f = constrain(constraint, f);
    
    // constraint.epsNativeSlice("paper/figs/rods-in-water-constraint.eps",
    // 			      Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
    // 			      Cartesian(0,ymax/2,zmax/2));
    //printf("Constraint has become a graph!\n");
   
    potential = water_prop.liquid_density*constraint
      + 100*water_prop.vapor_density*VectorXd::Ones(gd.NxNyNz);
    //potential = water_prop.liquid_density*VectorXd::Ones(gd.NxNyNz);
    potential = -water_prop.kT*potential.cwise().log();
    
    // {
    //    fflush(stdout);
    //    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    //    //density.epsNative1d("paper/figs/constrained-water-1D.eps",
    //    //			Cartesian(0,0,0), Cartesian(0,0,zmax),
    //    //			water_prop.liquid_density, 1, " ");
    
    //    //sleep(3);
    //  }
    
    Minimizer min = Precision(0, PreconditionedConjugateGradient(f, gd, water_prop.kT, 
								 &potential,
								 QuadraticLineMinimizer));
    
    printf("\nDistance between rods = %g bohr (%g nm)\n", distance, distance/nm);
    
    const int numiters = 200;
    for (int i=0;i<numiters && min.improve_energy(true);i++) {
      fflush(stdout);
      Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
      //density.epsNative1d("paper/figs/constrained-water-1D.eps",
      //			Cartesian(0,0,0), Cartesian(0,0,zmax),
      //			water_prop.liquid_density, 1, " ");
      
      //density.epsNativeSlice("paper/figs/rods-in-water.eps", 
      //			     Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
      //		     Cartesian(0,ymax/2,zmax/2));
      
      //sleep(3);
    }
    //min.print_info();

    // char *plotnamedens = (char *)malloc(1024);
    // sprintf(plotnamedens, "paper/figs/rods-density-1nm-%04.1f.dat", distance/nm);
    // density.epsNativeSlice(plotnamedens, Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
    // 			   Cartesian(0,ymax/2,zmax/2));
    // free(plotnamedens);

    double energy = min.energy()/width;
    printf("Energy is %.15g\n", energy);

    fprintf(o, "%g\t%.15g\n", distance/nm, energy);

    char *plotnameslice = (char *)malloc(1024);
    sprintf(plotnameslice, "paper/figs/rods-slice-1nm-%04.1f.dat", distance/nm);
    char *plotname = (char *)malloc(1024);
    sprintf(plotname, "paper/figs/rods-1nm-%04.1f.dat", distance/nm);
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    Grid energy_density(gd, f(water_prop.kT, gd, potential));
    Grid entropy(gd, S(water_prop.kT, potential));
    Grid Xassoc(gd, X(water_prop.kT, density));
    plot_grids_y_direction(plotnameslice, density, energy_density, entropy, Xassoc);
    plot_grids_yz_directions(plotname, density, 
			     energy_density, entropy, Xassoc);
    free(plotnameslice);
    free(plotname);

    // double N = 0;
    // {
    //   Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    //   for (int i=0;i<gd.NxNyNz;i++) N += density[i]*gd.dvolume;
    // }
    
    //N = N/width;
    //printf("N is %.15g\n", N);
  }
  fclose(o);
}
