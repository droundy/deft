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

static const double lj_pressure = 100000*3.3989316e-14; // 0.1 MPa in Hartree/bohr^3
static const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin
const double nm = 18.8972613;
const double angstrom = .1*nm;
const double kJpermol = 1/2625.49962; //kiloJoule per mole in hartree
// Here we set up the lattice.
double zmax = 2.5*nm;
double ymax = 2.5*nm;
double xmax = 2.5*nm;
double temperature; // temperature in Hartree
const double sigma = 3.57*angstrom; // in Bohr radii Xe in water
//from Dzublella, Swanson, and McCammon
const double epsilon = 1.0773*kJpermol; // in Hartree for Xe in water
double diameter = 0.714*nm;

double externalpotentialfunction(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  const double x = r.x();
  const double dist = sqrt(x*x+y*y+z*z);
  const double oodist6 = 1.0/uipow(dist/sigma, 6);
  const double pot = 4*epsilon*(oodist6*oodist6 - oodist6);
  const double max_pot = 5*temperature;
  if (pot < max_pot) return pot;
  return max_pot;
}

double notinwall(Cartesian r) {
  const double z = r.z();
  const double y = r.y();
  const double x = r.x();
  if (sqrt(sqr(z)+sqr(y)+sqr(x)) < diameter/2) {
      return 0; 
  }
  return 1;
}

void plot_grids_y_direction(const char *fname, const Grid &a, const Grid &b) {
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
    fprintf(out, "%g\t%g\t%g\t%g\t%g\n", here[0], here[1], here[2], ahere, bhere);
  }
  fclose(out);
}

int main(int argc, char *argv[]) {
  clock_t start_time = clock();
  if (argc == 2) {
    if (sscanf(argv[1], "%lg", &temperature) != 1) {
      printf("Got bad argument: %s\n", argv[1]);
      return 1;
    }
    temperature *= kB;

  } else {
    printf("Need temperature.\n");
    return 1;
  }

  printf("Diameter is %g bohr = %g nm\n", diameter, diameter/nm);
  const double padding = 1*nm;
  // xmax = ymax = zmax = diameter + 2*padding;

  Functional f_sphere = OfEffectivePotential(SaftFluid2(hughes_water_prop.lengthscale,
						hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
						hughes_water_prop.epsilon_dispersion,
						hughes_water_prop.lambda_dispersion,
						hughes_water_prop.length_scaling, 0));
  double n_sphere = pressure_to_density(f_sphere, temperature, lj_pressure,
				      0.001, 0.01);

  double mu_sphere = find_chemical_potential(f_sphere, temperature, n_sphere);

  f_sphere = OfEffectivePotential(SaftFluid2(hughes_water_prop.lengthscale,
				     hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
				     hughes_water_prop.epsilon_dispersion,
				     hughes_water_prop.lambda_dispersion,
				     hughes_water_prop.length_scaling, mu_sphere));

  Functional X_sphere = WaterX(hughes_water_prop.lengthscale,
                        hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
                        hughes_water_prop.epsilon_dispersion,
                        hughes_water_prop.lambda_dispersion,
                        hughes_water_prop.length_scaling, mu_sphere);
  
  Functional HB_sphere = HughesHB(hughes_water_prop.lengthscale,
                        hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
                        hughes_water_prop.epsilon_dispersion,
                        hughes_water_prop.lambda_dispersion,
                        hughes_water_prop.length_scaling, mu_sphere);

  Functional S_sphere = OfEffectivePotential(EntropySaftFluid2(hughes_water_prop.lengthscale,
                                                        hughes_water_prop.epsilonAB,
                                                        hughes_water_prop.kappaAB,
                                                        hughes_water_prop.epsilon_dispersion,
                                                        hughes_water_prop.lambda_dispersion,
                                                        hughes_water_prop.length_scaling));
  
  const double EperVolume_sphere = f_sphere(temperature, -temperature*log(n_sphere));
  const double EperNumber_sphere = EperVolume_sphere/n_sphere;
  const double SperNumber_sphere = S_sphere(temperature, -temperature*log(n_sphere))/n_sphere;
  const double EperCell_sphere = EperVolume_sphere*(zmax*ymax*xmax - (M_PI/6)*diameter*diameter*diameter);

  Lattice lat(Cartesian(xmax,0,0), Cartesian(0,ymax,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, 0.1);

  Grid sphere_potential(gd);
  Grid sphere_constraint(gd);
  sphere_constraint.Set(notinwall);

  f_sphere = OfEffectivePotential(SaftFluid2(hughes_water_prop.lengthscale,
				       hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
				       hughes_water_prop.epsilon_dispersion,
				       hughes_water_prop.lambda_dispersion,
				       hughes_water_prop.length_scaling, mu_sphere));
  f_sphere = constrain(sphere_constraint, f_sphere);

  sphere_potential = hughes_water_prop.liquid_density*sphere_constraint
      + 100*hughes_water_prop.vapor_density*VectorXd::Ones(gd.NxNyNz);
  sphere_potential = -temperature*sphere_potential.cwise().log();

  double energy_sphere;
    {
      const double surface_tension = 5e-5; // crude guess from memory...
      const double surfprecision = 1e-4*M_PI*diameter*diameter*surface_tension; // four digits precision
      const double bulkprecision = 1e-12*fabs(EperCell_sphere); // but there's a limit on our precision for small spheres
      const double precision = bulkprecision + surfprecision;
      Minimizer min = Precision(precision,
                                PreconditionedConjugateGradient(f_sphere, gd, temperature, 
                                                                &sphere_potential,
                                                                QuadraticLineMinimizer));

      printf("\nDiameter of sphere = %g bohr (%g nm)\n", diameter, diameter/nm);

      const int numiters = 200;
      for (int i=0;i<numiters && min.improve_energy(true);i++) {

        double peak = peak_memory()/1024.0/1024;
        double current = current_memory()/1024.0/1024;
        printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      }

      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      
      energy_sphere = min.energy();
      printf("Total energy is %.15g\n", energy_sphere);
      Grid sphere_density(gd, EffectivePotentialToDensity()(temperature, gd, sphere_potential));
      // Here we free the minimizer with its associated data structures.
    }

    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
    }

    double entropy_sphere = S_sphere.integral(temperature, sphere_potential);
    Grid sphere_density(gd, EffectivePotentialToDensity()(temperature, gd, sphere_potential));
    {
    printf("Number of water molecules is %g\n", sphere_density.integrate());
    printf("The bulk energy per cell should be %g\n", EperCell_sphere);
    printf("The bulk energy based on number should be %g\n", EperNumber_sphere*sphere_density.integrate());
    printf("The bulk entropy is %g/N\n", SperNumber_sphere);
    Functional otherS_sphere = EntropySaftFluid2(hughes_water_prop.lengthscale,
                                          hughes_water_prop.epsilonAB,
                                          hughes_water_prop.kappaAB,
                                          hughes_water_prop.epsilon_dispersion,
                                          hughes_water_prop.lambda_dispersion,
                                          hughes_water_prop.length_scaling);
    printf("The bulk entropy (haskell) = %g/N\n", otherS_sphere(temperature, n_sphere)/n_sphere);
    //printf("My entropy is %g when I would expect %g\n", entropy, entropy - SperNumber*sphere_density.integrate());
    double hentropy_sphere = otherS_sphere.integral(temperature, sphere_density);
    otherS_sphere.print_summary("   ", hentropy_sphere, "total entropy");
    printf("My haskell entropy is %g, when I would expect = %g, difference is %g\n", hentropy_sphere,
           otherS_sphere(temperature, n_sphere)*sphere_density.integrate()/n_sphere,
           hentropy_sphere - otherS_sphere(temperature, n_sphere)*sphere_density.integrate()/n_sphere);

    Grid zeroed_out_density(gd, sphere_density.cwise()*sphere_constraint); // this is zero inside the sphre!  :
    //char *plotname = (char *)malloc(1024);
    //sprintf(plotname, "papers/water-saft/figs/hughes-sphere-%04.2f.dat", diameter/nm);
    //plot_grids_y_direction(plotname, sphere_density, X_values, H_bonds);

    //free(plotname);

    //sphere_density.epsNativeSlice("papers/hughes-saft/figs/sphere.eps",
		//	   Cartesian(0,ymax,0), Cartesian(0,0,zmax),
		//	   Cartesian(0,ymax/2,zmax/2));
    
    double peak = peak_memory()/1024.0/1024;
    printf("Peak memory use is %g M\n", peak);

    double oldN = sphere_density.integrate();
    sphere_density = n_sphere*VectorXd::Ones(gd.NxNyNz);
    double hentropyb_sphere = otherS_sphere.integral(temperature, sphere_density);
    // printf("bulklike thingy has %g molecules\n", sphere_density.integrate());
    otherS_sphere.print_summary("   ", hentropyb_sphere, "bulk-like entropy");
    printf("entropy difference is %g\n", hentropy_sphere - hentropyb_sphere*oldN/sphere_density.integrate());
    }
  // }

    // FILE *o = fopen(datname, "w");
    // fprintf(o, "%g\t%.15g\n", diameter/nm, energy - EperCell);
    // fprintf(o, "%g\t%.15g\t%.15g\t%.15g\t%.15g\t%g\n", diameter/nm, energy - EperNumber*sphere_density.integrate(), energy - EperCell,
    //        temperature*(entropy - SperNumber*sphere_density.integrate()),
    //        temperature*(hentropy - otherS(temperature, n_sphere)*sphere_density.integrate()/n_sphere), broken_H_bonds);
    // fclose(o);

///////////////////Lennard-Jones stuff////////////////////////////////////

  char *datname = (char *)malloc(1024);
  sprintf(datname, "papers/water-saft/figs/hughes-lj-XE-%gK-energy-hs.dat", temperature/kB);

  Functional f = OfEffectivePotential(SaftFluid2(hughes_water_prop.lengthscale,
                                                hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
                                                hughes_water_prop.epsilon_dispersion,
                                                hughes_water_prop.lambda_dispersion,
                                                hughes_water_prop.length_scaling, 0));
  double n_lj = pressure_to_density(f, temperature, lj_pressure,
                                      0.001, 0.01);

  double mu_lj = find_chemical_potential(f, temperature, n_lj);

  f = OfEffectivePotential(SaftFluid2(hughes_water_prop.lengthscale,
                                     hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
                                     hughes_water_prop.epsilon_dispersion,
                                     hughes_water_prop.lambda_dispersion,
                                     hughes_water_prop.length_scaling, mu_lj));

  Functional S = OfEffectivePotential(EntropySaftFluid2(hughes_water_prop.lengthscale,
                                                        hughes_water_prop.epsilonAB,
                                                        hughes_water_prop.kappaAB,
                                                        hughes_water_prop.epsilon_dispersion,
                                                        hughes_water_prop.lambda_dispersion,
                                                        hughes_water_prop.length_scaling));
  
  const double EperVolume = f(temperature, -temperature*log(n_lj));
  const double EperNumber = EperVolume/n_lj;
  const double SperNumber = S(temperature, -temperature*log(n_lj))/n_lj;
  const double EperCell = EperVolume*(zmax*ymax*xmax - (M_PI/6)*sigma*sigma*sigma);
    
  Grid potential(gd);
  Grid externalpotential(gd);
  externalpotential.Set(externalpotentialfunction);
    
  f = OfEffectivePotential(SaftFluid2(hughes_water_prop.lengthscale,
                                     hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
                                     hughes_water_prop.epsilon_dispersion,
                                     hughes_water_prop.lambda_dispersion,
                                     hughes_water_prop.length_scaling, mu_lj) + ExternalPotential(externalpotential));

  Functional X = WaterX(hughes_water_prop.lengthscale,
                        hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
                        hughes_water_prop.epsilon_dispersion,
                        hughes_water_prop.lambda_dispersion,
                        hughes_water_prop.length_scaling, mu_lj);
  
  Functional HB = HughesHB(hughes_water_prop.lengthscale,
                           hughes_water_prop.epsilonAB, hughes_water_prop.kappaAB,
                           hughes_water_prop.epsilon_dispersion,
                           hughes_water_prop.lambda_dispersion,
                           hughes_water_prop.length_scaling, mu_lj);

  externalpotential.epsNativeSlice("papers/water-saft/figs/hughes-lj-potential-hs.eps",
                                   Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
                                   Cartesian(0,ymax/2,zmax/2));
  printf("Done outputting hughes-lj-potential-hs.eps\n");

  potential = sphere_density; // ???
    
  double energy;
  {
    const double surface_tension = 5e-5; // crude guess from memory...
    const double surfprecision = 1e-4*M_PI*sigma*sigma*surface_tension; // four digits precision
    const double bulkprecision = 1e-12*fabs(EperCell); // but there's a limit on our precision
    const double precision = bulkprecision + surfprecision;
    Minimizer min = Precision(precision,
                              PreconditionedConjugateGradient(f, gd, temperature, 
                                                              &potential,
                                                              QuadraticLineMinimizer));
      
      
    const int numiters = 200;
    for (int i=0;i<numiters && min.improve_energy(true);i++) {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      fflush(stdout);
      // char* name = new char[1000];
      // sprintf(name, "papers/water-saft/figs/lj-%s-%d-density-big.eps", argv[1], i);
      // Grid density(gd, EffectivePotentialToDensity()(temperature, gd, potential));
      // density.epsNativeSlice(name,
      //                        Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
      //                        Cartesian(0,ymax/2,zmax/2));
      // Grid gradient(gd, potential);
      // gradient *= 0;
      // f.integralgrad(temperature, potential, &gradient);
      // sprintf(name, "papers/water-saft/figs/lj-%s-%d-gradient-big.eps", argv[1], i);
      // gradient.epsNativeSlice("papers/water-saft/figs/lj-gradient-big.eps",
      //                         Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
      //                         Cartesian(0,ymax/2,zmax/2));
      // sprintf(name, "papers/water-saft/figs/lj-%s-%d-big.dat", argv[1], i);
      // plot_grids_y_direction(name, density, gradient);
    }
    {
      char* name = new char[1000];
      sprintf(name, "papers/water-saft/figs/hughes-lj-Xe-%gK-density-hs-pot.eps", temperature/kB);
      Grid density(gd, EffectivePotentialToDensity()(temperature, gd, potential));
      density.epsNativeSlice(name,
                             Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
                             Cartesian(0,ymax/2,zmax/2));
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

  Grid gradient(gd, potential);
  gradient *= 0;
  f.integralgrad(temperature, potential, &gradient);
  gradient.epsNativeSlice("papers/water-saft/figs/hughes-lj-gradient-hs-pot.eps",
                          Cartesian(0,ymax,0), Cartesian(0,0,zmax), 
                          Cartesian(0,ymax/2,zmax/2));

  double entropy = S.integral(temperature, potential);
  Grid density(gd, EffectivePotentialToDensity()(temperature, gd, potential));
  // Grid zeroed_out_density(gd, density.cwise()*sphere_constraint); // this is zero inside the sphere!
  Grid X_values(gd, X(temperature, gd, density));
  //Grid H_bonds_grid(gd, zeroed_out_density.cwise()*(4*(VectorXd::Ones(gd.NxNyNz)-X_values)));
  //const double broken_H_bonds = (HB(temperature, n_lj)/n_lj)*zeroed_out_density.integrate() - H_bonds_grid.integrate();
  //printf("Number of water molecules is %g\n", density.integrate());
  printf("The bulk energy per cell should be %g\n", EperCell);
  printf("The bulk energy based on number should be %g\n", EperNumber*density.integrate());
  printf("The bulk entropy is %g/N\n", SperNumber);
  Functional otherS = EntropySaftFluid2(hughes_water_prop.lengthscale,
                                        hughes_water_prop.epsilonAB,
                                        hughes_water_prop.kappaAB,
                                        hughes_water_prop.epsilon_dispersion,
                                        hughes_water_prop.lambda_dispersion,
                                        hughes_water_prop.length_scaling);
  printf("The bulk entropy (haskell) = %g/N\n", otherS(temperature, n_lj)/n_lj);
  //printf("My entropy is %g when I would expect %g\n", entropy, entropy - SperNumber*density.integrate());
  double hentropy = otherS.integral(temperature, density);
  otherS.print_summary("   ", hentropy, "total entropy");
  printf("My haskell entropy is %g, when I would expect = %g, difference is %g\n", hentropy,
         otherS(temperature, n_lj)*density.integrate()/n_lj,
         hentropy - otherS(temperature, n_lj)*density.integrate()/n_lj);

  FILE *o = fopen(datname, "w");
  fprintf(o, "%g\t%.15g\t%.15g\t%.15g\n", temperature/kB, energy - EperNumber*density.integrate(),
          temperature*(entropy - SperNumber*density.integrate()),
          temperature*(hentropy - otherS(temperature, n_lj)*density.integrate()/n_lj));
  fclose(o);

  char *plotname = (char *)malloc(1024);
  sprintf(plotname, "papers/water-saft/figs/hughes-lj-Xe-%gK-hs-pot.dat",  temperature/kB);
  plot_grids_y_direction(plotname, density, gradient);

  free(plotname);

  double peak = peak_memory()/1024.0/1024;
  printf("Peak memory use is %g M\n", peak);

  double oldN = density.integrate();
  density = n_lj*VectorXd::Ones(gd.NxNyNz);;
  double hentropyb = otherS.integral(temperature, density);
  printf("bulklike thingy has %g molecules\n", density.integrate());
  otherS.print_summary("   ", hentropyb, "bulk-like entropy");
  printf("entropy difference is %g\n", hentropy - hentropyb*oldN/density.integrate());

  clock_t end_time = clock();
  double seconds = (end_time - start_time)/double(CLOCKS_PER_SEC);
  double hours = seconds/60/60;
  printf("Entire calculation took %.0f hours %.0f minutes\n", hours, 60*(hours-floor(hours)));
}
