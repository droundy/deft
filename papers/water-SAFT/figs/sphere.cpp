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

const double nm = 18.8972613;
// Here we set up the lattice.
double zmax = 2.5*nm;
double ymax = 2.5*nm;
double xmax = 2.5*nm;
double diameter = 1*nm;
bool using_default_diameter = true;

double notinwall(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  const double x = r.x();
  if (sqrt(sqr(z)+sqr(y)+sqr(x)) < diameter/2) {
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
  const int x = 0;
  const int z = 0;
  for (int y=0; y<gd.Ny/2; y++) {
    Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
    double ahere = a(x,y,z);
    fprintf(out, "%g\t%g\t%g\t%g\n", here[0], here[1], here[2], ahere);
  }
  fclose(out);
}

int main(int argc, char *argv[]) {
  clock_t start_time = clock();
  if (argc > 1) {
    if (sscanf(argv[1], "%lg", &diameter) != 1) {
      printf("Got bad argument: %s\n", argv[1]);
      return 1;
    }
    diameter *= nm;
    using_default_diameter = false;
  }
  printf("Diameter is %g bohr = %g nm\n", diameter, diameter/nm);
  const double padding = 1*nm;
  xmax = ymax = zmax = diameter + 2*padding;

  char *datname = (char *)malloc(1024);
  sprintf(datname, "papers/water-SAFT/figs/sphere-%04.2fnm-energy.dat", diameter/nm);
  
  Functional f = OfEffectivePotential(SaftFluid2(water_prop.lengthscale,
						water_prop.epsilonAB, water_prop.kappaAB,
						water_prop.epsilon_dispersion,
						water_prop.lambda_dispersion,
						water_prop.length_scaling, 0));
  double n_1atm = pressure_to_density(f, water_prop.kT, atmospheric_pressure,
				      0.001, 0.01);

  double mu_satp = find_chemical_potential(f, water_prop.kT, n_1atm);

  f = OfEffectivePotential(SaftFluid2(water_prop.lengthscale,
				     water_prop.epsilonAB, water_prop.kappaAB,
				     water_prop.epsilon_dispersion,
				     water_prop.lambda_dispersion,
				     water_prop.length_scaling, mu_satp));
  
  const double EperVolume = f(water_prop.kT, -water_prop.kT*log(n_1atm));
  const double EperNumber = EperVolume/n_1atm;
  const double EperCell = EperVolume*(zmax*ymax*xmax - (M_PI/6)*diameter*diameter*diameter);

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
  
  //for (diameter=0*nm; diameter<3.0*nm; diameter+= .1*nm) {
    Lattice lat(Cartesian(xmax,0,0), Cartesian(0,ymax,0), Cartesian(0,0,zmax));
    GridDescription gd(lat, 0.2);
    
    Grid potential(gd);
    Grid constraint(gd);
    constraint.Set(notinwall);
    
    f = OfEffectivePotential(SaftFluid2(water_prop.lengthscale,
				       water_prop.epsilonAB, water_prop.kappaAB,
				       water_prop.epsilon_dispersion,
				       water_prop.lambda_dispersion,
				       water_prop.length_scaling, mu_satp));
    f = constrain(constraint, f);
    //constraint.epsNativeSlice("papers/water-SAFT/figs/sphere-constraint.eps",
    // 			      Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
    // 			      Cartesian(0,ymax/2,zmax/2));
    //printf("Constraint has become a graph!\n");
   
    potential = water_prop.liquid_density*constraint
      + 100*water_prop.vapor_density*VectorXd::Ones(gd.NxNyNz);
    //potential = water_prop.liquid_density*VectorXd::Ones(gd.NxNyNz);
    potential = -water_prop.kT*potential.cwise().log();
    
    double energy;
    {
      const double surface_tension = 5e-5; // crude guess from memory...
      const double surfprecision = 1e-4*M_PI*diameter*diameter*surface_tension; // four digits precision
      const double bulkprecision = 1e-12*fabs(EperCell); // but there's a limit on our precision for small spheres
      const double precision = bulkprecision + surfprecision;
      Minimizer min = Precision(precision,
                                PreconditionedConjugateGradient(f, gd, water_prop.kT, 
                                                                &potential,
                                                                QuadraticLineMinimizer));
      
      printf("\nDiameter of sphere = %g bohr (%g nm)\n", diameter, diameter/nm);
      
      const int numiters = 200;
      for (int i=0;i<numiters && min.improve_energy(true);i++) {
        //fflush(stdout);
        //Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
        
        //density.epsNativeSlice("papers/water-SAFT/figs/sphere.eps", 
        //			     Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
        //			     Cartesian(0,ymax/2,zmax/2));
        
        //sleep(3);

        double peak = peak_memory()/1024.0/1024;
        double current = current_memory()/1024.0/1024;
        printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      }
      
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      
      energy = min.energy();
      printf("Total energy is %.15g\n", energy);
      // Here we free the minimizer with its associated data structures.
    }

    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
    }

    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    printf("Number of water molecules is %g\n", density.integrate());
    printf("The bulk energy per cell should be %g\n", EperCell);
    printf("The bulk energy based on number should be %g\n", EperNumber*density.integrate());

    FILE *o = fopen(datname, "w");
    //fprintf(o, "%g\t%.15g\n", diameter/nm, energy - EperCell);
    fprintf(o, "%g\t%.15g\t%.15g\n", diameter/nm, energy - EperNumber*density.integrate(), energy - EperCell);
    fclose(o);

    char *plotname = (char *)malloc(1024);

    sprintf(plotname, "papers/water-SAFT/figs/sphere-%04.2f.dat", diameter/nm);
    plot_grids_y_direction(plotname, density);

    free(plotname);

    //density.epsNativeSlice("papers/water-SAFT/figs/sphere.eps",
		//	   Cartesian(0,ymax,0), Cartesian(0,0,zmax),
		//	   Cartesian(0,ymax/2,zmax/2));
    
    double peak = peak_memory()/1024.0/1024;
    printf("Peak memory use is %g M\n", peak);
  
  // }
  clock_t end_time = clock();
  double seconds = (end_time - start_time)/double(CLOCKS_PER_SEC);
  double hours = seconds/60/60;
  printf("Entire calculation took %.0f hours %.0f minutes\n", hours, 60*(hours-floor(hours)));
}
