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
#include "Testables.h"
#include "equation-of-state.h"
#include "LineMinimizer.h"

#include "utilities.h"

double zmax = 150;
const double width = 0.0001;
double cavitysize = 60;

double notinwall(Cartesian r) {
  const double z = r.z();
  if (fabs(z) > (zmax-cavitysize)/2) {
    return 1;
  } else {
    return 0;
  }
}

void plot_grids_z_direction(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = a.description();
  const int x = gd.Nx/2;
  const int y = gd.Ny/2;
  for (int z=0; z<gd.Nz; z++) {
    Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
    double ahere = a(x,y,z);
    double bhere = b(x,y,z);
    double chere = c(x,y,z);
    double dhere = d(x,y,z);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[0], here[1], here[2], ahere, bhere, chere, dhere);
  }
  fclose(out);
}

int main(int, char **) {
  FILE *o = fopen("paper/figs/constrained-water-1D.dat", "w");

  Functional f = SaftFluid(water_prop.lengthscale,
			   water_prop.epsilonAB, water_prop.kappaAB,
			   water_prop.epsilon_dispersion,
			   water_prop.lambda_dispersion,
			   water_prop.length_scaling, 0);
  double mu_satp = find_chemical_potential(f, water_prop.kT,
                                           water_prop.liquid_density);
  f = SaftFluid(water_prop.lengthscale,
		water_prop.epsilonAB, water_prop.kappaAB,
		water_prop.epsilon_dispersion,
		water_prop.lambda_dispersion,
		water_prop.length_scaling, mu_satp);
    
  Functional X = Xassociation(water_prop.lengthscale, water_prop.epsilonAB, 
  			    water_prop.kappaAB, water_prop.epsilon_dispersion,
  			    water_prop.lambda_dispersion,
  			    water_prop.length_scaling);
  
  Functional S = SaftEntropy(water_prop.lengthscale, water_prop.epsilonAB, 
				   water_prop.kappaAB, water_prop.epsilon_dispersion,
				   water_prop.lambda_dispersion,
				   water_prop.length_scaling);
  
  for (cavitysize = 5; cavitysize<zmax; cavitysize+=5) {

    //printf("Cavity size is %g bohr\n", cavitysize);
    
    Lattice lat(Cartesian(width,0,0), Cartesian(0,width,0), Cartesian(0,0,zmax));
    GridDescription gd(lat, 0.2);
   
    Grid potential(gd);
    Grid constraint(gd);
    constraint.Set(notinwall);
    //f = constrain(constraint, f);
    
    potential = water_prop.liquid_density*constraint
      + 100*water_prop.vapor_density*VectorXd::Ones(gd.NxNyNz);
    //potential = water_prop.liquid_density*VectorXd::Ones(gd.NxNyNz);
    potential = -water_prop.kT*potential.cwise().log();
    
    Minimizer min = Precision(0, 
			      PreconditionedConjugateGradient(f, gd, water_prop.kT, 
							      &potential,
							      QuadraticLineMinimizer));
    printf("Cavity size is %g bohr\n", cavitysize);

    const int numiters = 100;
    for (int i=0;i<numiters && min.improve_energy(true);i++) {
      fflush(stdout);
      Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
      // density.epsNative1d("paper/figs/constrained-water-1D.eps",
      // 			  Cartesian(0,0,0), Cartesian(0,0,zmax),
      // 			  water_prop.liquid_density, 1, " ");
      //sleep(1);
    }
    //min.print_info();
    
    double energy = min.energy()/width/width;
    //printf("Energy is %.15g\n", energy);

    fprintf(o, "\t%g\t%.15g\n", cavitysize, energy);

    char *plotname = (char *)malloc(1024);
    sprintf(plotname, "paper/figs/cavitysize-%03g.dat", cavitysize);
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    Grid energy_density(gd, f(water_prop.kT, gd, potential));
    Grid entropy(gd, S(water_prop.kT, potential));
    Grid Xassoc(gd, X(water_prop.kT, density));

    // double radius = water_prop.lengthscale;
    // Functional R(radius, "R");
    // Functional n2 = ShellConvolve(radius);
    // Functional n0 = n2/(4*M_PI*sqr(R));
    // Functional delta = DeltaSAFT(water_prop.lengthscale, water_prop.epsilonAB,
    // 				 water_prop.kappaAB,
    // 				 water_prop.epsilon_dispersion,
    // 				 water_prop.lambda_dispersion,
    // 				 water_prop.length_scaling);
    // Functional zeta = getzeta(radius);
    // Functional fourn0zetadelta = 4* n0 * zeta*delta;
    // Grid foobar(gd, fourn0zetadelta(water_prop.kT, density));
    plot_grids_z_direction(plotname, density, energy_density, entropy, Xassoc);
    //plot_grids_z_direction(plotname, density, energy_density, foobar, Xassoc);
    free(plotname);
    
    //double N = 0;
    //{
    // Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    // for (int i=0;i<gd.NxNyNz;i++) N += density[i]*gd.dvolume;
    //}
    //N = N/width/width;
    //printf("N is %.15g\n", N);
    
    //Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    density.epsNative1d("paper/figs/constrained-water-1D.eps", 
    			Cartesian(0,0,0), Cartesian(0,0,zmax), 
    			water_prop.liquid_density, 1, " ");
    //potential.epsNative1d("hard-wall-potential.eps", Cartesian(0,0,0), Cartesian(0,0,zmax), 1, 1);
   
  }
  fclose(o);
}
