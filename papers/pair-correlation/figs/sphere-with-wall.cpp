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
#include <math.h>
#include "OptimizedFunctionals.h"
#include "equation-of-state.h"
#include "LineMinimizer.h"
#include "ContactDensity.h"
#include "utilities.h"
#include "handymath.h"
#include "errno.h"
#include "sys/stat.h" // for mkdir

int count =0;

// Maximum and spacing values for plotting

// Here we set up the lattice.
static double width = 20;
const double dw = 0.001;
const double dx = 0.1;
const double spacing = 3; // space on each side


double notinwall_or_sphere(Cartesian r) {
  const double x = r.x();
  const double y = r.y();
  const double z = r.z();
  if (z < spacing || width-z < spacing) {
    return 0;
  }
  const double dis = sqrt((z-spacing)*(z-spacing)+y*y*+x*x);
  if (dis<2.0) {
    return 0;
  }
  return 1;
}

void pair_plot(const char *fname, const char *fname_debug, const Grid &density) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  FILE *out_debug = fopen(fname_debug, "w");
  if (!out_debug) {
    fprintf(stderr, "Unable to create file %s!\n", fname_debug);
    // don't just abort?
    return;
  }
  printf("In function one\n");
  fflush(stdout);
  const GridDescription gd = density.description();
  const int x = 0;
  printf("In function two\n");
  fflush(stdout);
  for (int y=0; y<gd.Nz/2; y++) {
    for (int z= int(spacing/dx+0.5); z<gd.Nz/2; z++) {
      printf("In function three\n");
      fflush(stdout);
      Cartesian here = gd.fineLat.toCartesian(Relative(x,y,z));
      printf("In function four\n");
      fflush(stdout);
      double pair_correlation = density(x,y,z)/density(x,y,gd.Nz-z);/*this is density
      on the sphere side divided by the density the same distance from
      the wall on the other side (where there is no sphere)*/
      printf("In function five\n");
      fflush(stdout);
      fprintf(out_debug, "%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[2], here[1], density(x,y,z), density(x,y,gd.Nz-z),
              density(x,y,-z), density(x,y,gd.Nz/2-z), pair_correlation);
      printf("In function six\n");
      fflush(stdout);
      fprintf(out, "%g\t", pair_correlation);
    }
    fprintf(out,"\n");
  }
  fclose(out_debug);
  fclose(out);
}

void path_plot(const char *fname, const Grid &density, const Grid &constraint) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
    return;
  }
  const GridDescription gd = density.description();
  const double z0 = spacing + 0.05;
  double radius_path = 2.005;
  int num = 100;
  double g2_path;
  double x_path;
  for (int i=0; i<int((10.0-radius_path)/dx+0.5); i++){
    x_path = i*dx;
    Cartesian r = Cartesian(10.0-x_path,0,z0);
    Cartesian r_opp_side = Cartesian(10.0-x_path,0,gd.Nz*dx-z0);
    g2_path = density(r)*constraint(r_opp_side)/density(r_opp_side)/constraint(r);
    fprintf(out,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x_path,g2_path,density(r),z0,constraint(r_opp_side),density(r_opp_side),constraint(r));
  }
  for (int i=0; i<num ;i++){
    double theta = i*M_PI/num/2.0;
    x_path = i*M_PI/num/2.0*radius_path + 8.0;
    Cartesian r = Cartesian(radius_path*cos(theta),0,z0+radius_path*sin(theta));
    Cartesian r_opp_side = Cartesian(radius_path*cos(theta),0,gd.Nz*dx-z0-radius_path*sin(theta));
    g2_path = density(r)*constraint(r_opp_side)/density(r_opp_side)/constraint(r);
    fprintf(out,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x_path,g2_path,density(r),z0,constraint(r_opp_side),density(r_opp_side),constraint(r));
  }
  for (int i=0; i<int(8.0/dx+0.5);i++){
    double r1z = i*dx;
    x_path = i*dx + radius_path*M_PI/2.0 + 10.0-radius_path;
    Cartesian r = Cartesian(0,0,r1z+radius_path+z0);
    Cartesian r_opp_side = Cartesian(0,0,gd.Nz*dx-r1z-radius_path-z0);
    g2_path = density(r)*constraint(r_opp_side)/density(r_opp_side)/constraint(r);
    fprintf(out,"%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x_path,g2_path,density(r),z0,constraint(r_opp_side),density(r_opp_side),constraint(r));
  }
  fclose(out);
}

Functional WB = HardSpheresNoTensor2(1.0);

int main(int, char **) {
  for (double eta = 0.3; eta < 0.35; eta += 0.1) {
    // Generates a data file for the pair distribution function, for filling fraction eta
    // and distance of first sphere from wall of z0. Data saved in a table such that the
    // columns are x values and rows are z1 values.
    printf("Now starting sphere_with_wall with eta = %g\n",eta);
    Lattice lat(Cartesian(width,0,0), Cartesian(0,width,0), Cartesian(0,0,width+2*spacing));
    printf("Hello???\n");
    GridDescription gd(lat, dx); //the resolution here dramatically affects our memory use
    printf("gd.Nz/2 = %d\n",gd.Nz/2);
    fflush(stdout);

    Functional f = OfEffectivePotential(WB + IdealGas());
    double mu = find_chemical_potential(f, 1, eta/(4*M_PI/3));
    f = OfEffectivePotential(WB + IdealGas()
                             + ChemicalPotential(mu));

    Grid potential(gd);
    Grid constraint(gd);

    constraint.Set(*notinwall_or_sphere);
    f = constrain(constraint, f);

    potential = (eta*constraint + 1e-4*eta*VectorXd::Ones(gd.NxNyNz))/(4*M_PI/3);
    potential = -potential.cwise().log();

    const double approx_energy = (WB + IdealGas() + ChemicalPotential(mu))(1, eta/(4*M_PI/3))*dw*dw*width;
    const double precision = fabs(approx_energy*1e-4);
    //printf("Minimizing to %g absolute precision...\n", precision);
    Minimizer min = Precision(precision,
                              PreconditionedConjugateGradient(f, gd, 1,
                                                              &potential,
                                                              QuadraticLineMinimizer));
    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      fflush(stdout);
    }

    for (int i=0;min.improve_energy(true) && i<100;i++) {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      fflush(stdout);
    }
    printf("HELLO!!one\n");
    fflush(stdout);
    Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
    printf("HELLO!!two\n");
    fflush(stdout);

    char *plotname = new char[1024];
    sprintf(plotname, "papers/pair-correlation/figs/wallsWB-with-sphere-%04.2f.dat", eta);
    printf("HELLO!!three\n");
    fflush(stdout);
    char *plotname_debug = new char[1024];
    sprintf(plotname_debug, "papers/pair-correlation/figs/wallsWB-with-sphere-debug-%04.2f.dat", eta);
    printf("HELLO!!four\n");
    fflush(stdout);
    pair_plot(plotname, plotname_debug, density);
    printf("HELLO!!five\n");
    fflush(stdout);
    delete[] plotname;
    delete[] plotname_debug;
    printf("HELLO!!six\n");
    fflush(stdout);
    char *plotname_path = new char[1024];
    sprintf(plotname_path, "papers/pair-correlation/figs/wallsWB-with-sphere-path-%04.2f.dat", eta);
    printf("HELLO!!seven\n");
    fflush(stdout);
    path_plot(plotname_path, density, constraint);
    printf("HELLO!!eight\n");
    fflush(stdout);
    delete[] plotname_path;
    printf("HELLO!!nine\n");
    fflush(stdout);
    
    {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      fflush(stdout);
    }
  }

  // Just create this file so make knows we have run.
  if (!fopen("papers/pair-correlation/figs/walls_sphere.dat", "w")) {
    printf("Error creating walls.dat!\n");
    return 1;
  }
  return 1;
}
