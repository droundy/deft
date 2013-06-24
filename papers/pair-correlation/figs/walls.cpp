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
double radial_distribution(double gsigma, double r);
double py_rdf (double eta, double r);
double mc (double local_eta, double r, double r_step, double g[]);


// Maximum and spacing values for plotting
const double zmax = 20;
const double xmax = 10;
const double dx = 0.1;

double mc_r_step;
double num_eta = 13.0;
double eta_step = 0.65/num_eta;
double * g = new double[1300*(int(num_eta)+1)];

// The functions for different ways of computing the pair distribution function.
double pairdist_this_work(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) {
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  return (radial_distribution(gsigma(r0), r) + radial_distribution(gsigma(r1), r))/2;
}
double pairdist_this_work_mc(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) {
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  const double gs0 = gsigma(r0);
  const double gs1 = gsigma(r1);
  const double c0 = -pow(fabs(sqrt(3.0)*sqrt(27.0*gs0*gs0*gs0*gs0 - 2.0*gs0*gs0*gs0)-9.0*gs0*gs0), 1.0/3.0);
  const double eta0 = c0/(pow(6.0, 2.0/3.0)*gs0) + 1.0/(pow(6.0, 1.0/3.0)*c0) + 1.0;
  const double c1 = -pow(fabs(sqrt(3.0)*sqrt(27.0*gs1*gs1*gs1*gs1 - 2.0*gs1*gs1*gs1)-9.0*gs1*gs1), 1.0/3.0);
  const double eta1 = c1/(pow(6.0, 2.0/3.0)*gs1) + 1.0/(pow(6.0, 1.0/3.0)*c1) + 1.0;
  //printf("gsigma0: %g, gsigma1: %g, eta0: %g, eta1: %g, c0: %g\n", gs0, gs1, eta0, eta1, c0);
  //fflush(stdout);
  return (mc(eta0, r, mc_r_step, g) + mc(eta1, r, mc_r_step, g))/2;
}

double pairdist_gross(const Grid &gsigma, const Grid &n, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) {
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  const double eta = 4.0/3*M_PI*1*1*1*(n(r0) + n(r1))/2;
  return mc(eta, r,mc_r_step,g);
}
double pairdist_fischer(const Grid &gsigma, const Grid &n, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) {
  // This implements the pair distribution function of Fischer and
  // Methfessel from the 1980 paper.  The py_rdf below should be the
  // true radial distribution function for a homogeneous hard-sphere
  // fluid with packing fraction eta.
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  const double eta = n3(Cartesian(0.5*(r0+r1)));
  return mc(eta, r, mc_r_step, g);
}

const char *fun[] = {
  "this-work",
  "this-work-mc",
  "gross",
  "fischer"
};
double (*pairdists[])(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) = {
  pairdist_this_work,
  pairdist_this_work_mc,
  pairdist_gross,
  pairdist_fischer
};
const int numplots = sizeof fun/sizeof fun[0];


// Here we set up the lattice.
static double width = 30;
const double dw = 0.001;
const double spacing = 3; // space on each side

double radial_distribution(double gsigma, double r) {
  // Constants determined by fit to monte-carlo data by plot-ghs.py
  const double a0 = 6.203,
               a1 = .4154,
               a2 = 6.449,
               a3 = 2.061,
               a4 = .3287,
               a5 = 4.555,
               a6 = 3.312;
  // const double a0 = 10.117,
  //              a1 = 0.309,
  //              a2 = 7.378,
  //              a3 = 1.569,
  //              a4 = 0.098,
  //              a5 = 6.053,
  //              a6 = 2.858;
  double d = r/2 - 1;
  if (d < 0)
    return 0;
  return 1 + (gsigma-1)*exp(-a0*d) + a1*(gsigma-1)*sin(a2*d)*exp(-a3*d)
           - a4*(gsigma-1)*(gsigma-1)*sin(a5*d)*exp(-a6*d);
}

void read_mc() {
  char *fname = new char[1024];
  sprintf(fname, "papers/pair-correlation/figs/gr-0.10.dat");
  FILE *in = fopen(fname, "r");
  if (!in) {
    fprintf(stderr, "Unable to open file %s!!!!\n", fname);
    delete[] fname;
    fclose(in);
    exit(1);
  }
  double num1;
  double num2;
  if (fscanf(in, " %lg %*g %lg", &num1, &num2) != 2) {
    printf("Error reading file...\n");
    delete[] fname;
    fclose(in);
    exit(1);
  }
  mc_r_step = 2*(num2 - num1);
  delete[] fname;
  fclose(in);
  for (int i=0; i<1300; i++) {
    g[i] = 1.0;
  }
  for (double eta = eta_step; eta < 0.5001; eta += eta_step) {
    char *fname = new char[1024];
    sprintf(fname, "papers/pair-correlation/figs/gr-%2.2f.dat", eta);
    FILE *in = fopen(fname, "r");
    if (!in) {
      fprintf(stderr, "Unable to open file %s!!!!\n", fname);
      delete[] fname;
      exit(1);
    }
    int i=0;
    while (fscanf(in, " %*g %lg", &num1) == 1) {
      g[int(floor(1300*(eta/eta_step)))+i] = num1/eta;
      i++;
    }
    if (i!=1300) {
      printf("There  are not 1300 lines in papers/pair-correlation/figs/gr-%2.2f.dat which is a problem in walls.cpp\n", eta);
      exit(1);
    }
    delete[] fname;
    fclose(in);
  }
}

double mc (double local_eta, double r, double r_step, double g[]) {
  int low_eta_i = floor(local_eta/eta_step);
  int high_eta_i = low_eta_i + 1;
  double g_low_eta;
  double g_high_eta;
  double fac;
  int j=0;
  double r_floor;
  if (local_eta > num_eta*eta_step) {
    //fprintf(stderr, "Error. Trying to interpolate to eta = %g, but our max is %g\n", local_eta, num_eta*eta_step);
    return NAN;
  }
  if (r < 2) {
    return 0;
  }
  if (r < 2 +(r_step/2.0)) {
    fac = (r-2.0)/(0.5*r_step);
    g_low_eta = (1-fac)*g[low_eta_i*1300]+fac*g[low_eta_i*1300 + 1];
    g_high_eta = (1-fac)*g[high_eta_i*1300]+fac*g[high_eta_i*1300 + 1];
  } else if (r > 14.985) {
    return 1.0;
  } else {
    j = floor((r-2.0+0.5*r_step)/r_step);
    r_floor = 2+(j-0.5)*r_step;
    fac = (r-r_floor)/r_step;
    g_low_eta = (1-fac)*g[low_eta_i*1300+j] + fac*g[low_eta_i*1300+j+1];
    g_high_eta = (1-fac)*g[high_eta_i*1300+j] + fac*g[high_eta_i*1300+j+1];
  }
  //printf("low_eta_i = %d and high_eta_i = %d\n", low_eta_i*1300+j, high_eta_i*1300+j);
  //fflush(stdout);
  fac = (local_eta-(low_eta_i)*eta_step)/eta_step;
  double ghere = (1-fac)*g_low_eta + fac*g_high_eta;
  return (1-fac)*g_low_eta + fac*g_high_eta;
}

double notinwall(Cartesian r) {
  const double z = r.z();
  if (fabs(z) > spacing) {
    return 1;
  }
  return 0;
}

static void took(const char *name) {
//   assert(name); // so it'll count as being used...
//   static clock_t last_time = clock();
//   clock_t t = clock();
//   double peak = peak_memory()/1024.0/1024;
//   double seconds = (t-last_time)/double(CLOCKS_PER_SEC);
//   if (seconds > 120) {
//     printf("\t\t%s took %.0f minutes and %g M memory\n", name, seconds/60, peak);
//   } else {
//     printf("\t\t%s took %g seconds and %g M memory\n", name, seconds, peak);
//   }
//   fflush(stdout);
//   last_time = t;
}

Functional WB = HardSpheresNoTensor2(1.0);

void z_plot(const char *fname, const Grid &a, const Grid &b, const Grid &c) {
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
    double bhere = b(x,y,z);
    double chere = c(x,y,z);
    fprintf(out, "%g\t%g\t%g\t%g\n", here[2], ahere, bhere, chere);
  }
  fclose(out);
}


void run_walls(double eta, const char *name, Functional fhs) {
  // Generates a data file for the pair distribution function, for filling fraction eta
  // and distance of first sphere from wall of z0. Data saved in a table such that the
  // columns are x values and rows are z1 values.
  printf("Now starting run_walls with eta = %g name = %s\n",
         eta, name);

  Functional f = OfEffectivePotential(fhs + IdealGas());
  double mu = find_chemical_potential(f, 1, eta/(4*M_PI/3));
  f = OfEffectivePotential(fhs + IdealGas()
                           + ChemicalPotential(mu));

  Lattice lat(Cartesian(dw,0,0), Cartesian(0,dw,0), Cartesian(0,0,width+2*spacing));
  GridDescription gd(lat, 0.01);

  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
  f = constrain(constraint, f);

  potential = (eta*constraint + 1e-4*eta*VectorXd::Ones(gd.NxNyNz))/(4*M_PI/3);
  potential = -potential.cwise().log();

  const double approx_energy = (fhs + IdealGas() + ChemicalPotential(mu))(1, eta/(4*M_PI/3))*dw*dw*width;
  const double precision = fabs(approx_energy*1e-4);
  //printf("Minimizing to %g absolute precision...\n", precision);
  Minimizer min = Precision(precision,
                            PreconditionedConjugateGradient(f, gd, 1,
                                                            &potential,
                                                            QuadraticLineMinimizer));
  for (int i=0;min.improve_energy(false) && i<100;i++) {
  }
  took("Doing the minimization");

  Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
  //printf("# per area is %g at filling fraction %g\n", density.sum()*gd.dvolume/dw/dw, eta);

  char *plotname = new char[1024];
  Grid gsigma(gd, gSigmaA(1.0)(1, gd, density));
  Grid nA(gd, ShellConvolve(2)(1, density)/(4*M_PI*4));
  Grid n3(gd, StepConvolve(1)(1, density));

  sprintf(plotname, "papers/pair-correlation/figs/walls%s-%04.2f.dat", name, eta);
  z_plot(plotname, density, gsigma, nA);

  // Create the walls directory if it doesn't exist.
  if (mkdir("papers/pair-correlation/figs/walls", 0777) != 0 && errno != EEXIST) {
    // We failed to create the directory, and it doesn't exist.
    printf("Failed to create papers/pair-correlation/figs/walls: %s",
           strerror(errno));
    exit(1); // fail immediately with error code
  }
  // here you choose the values of z0 to use
  //dx is set at beggining of file
  for (double z0 = 3.05; z0 < 13; z0 += dx) {
    // For each z0, we now pick one of our methods for computing the
    // pair distribution function:
    for (int version = 0; version < numplots; version++) {
      sprintf(plotname,
              "papers/pair-correlation/figs/walls/walls%s-%s-pair-%04.2f-%1.2f.dat",
              name, fun[version], eta, z0-3);
      FILE *out = fopen(plotname, "w");
      if (!out) {
        fprintf(stderr, "Unable to create file %s!\n", plotname);
        return;
      }
      // the +1 for z0 and z1 are to shift the plot over, so that a sphere touching the wall
      // is at z = 0, to match with the monte carlo data
      const Cartesian r0(0,0,z0);
      for (double x = 0; x < xmax - dx/2; x += dx) {
        for (double z1 = 3; z1 < zmax + 3 - dx/2; z1 += dx) {
          const Cartesian r1(x,0,z1);
          double g2 = pairdists[version](gsigma, density, nA, n3, r0, r1);
          fprintf(out, "%g\t", g2);
        }
        fprintf(out, "\n");
      }
      fclose(out);
      took(plotname);
    }
  }
  delete[] plotname;
  took("Dumping the pair dist plots");
  //This is the begginning of the integral to get a1.  It takes way to long (finished about an eighth of it when I left
  //it running over night.  So I'm wondering - Can this be fixed? Should we put it in a different file?
  took("Making 2d plots");
  //printf("Starting the a1 integrals now!!\n");
  for (int version = 0; version < numplots; version++) {
    for (double delta_r = 2.0; delta_r <= 4.0; delta_r += 0.1){
      const double dv = 0.01;
      char *plotname_a = new char [1024];
      sprintf(plotname_a, "papers/pair-correlation/figs/walls/walls_da%s-%s-%04.2f-%04.2f.dat", name, fun[version], eta, delta_r);
      FILE *out = fopen(plotname_a,"w");
      if (!out) {
        fprintf(stderr, "Unable to create file %s!\n", plotname_a);
        return;
      }
      delete[] plotname_a;
      for (double z0 = 2; z0 < 13; z0 += dx) {
        double da_dz = 0;
        const Cartesian r0(0,0,z0);
        const double dtheta = M_PI/ceil(delta_r/dv*M_PI);
        for (double theta = dtheta/2; theta <= M_PI; theta += dtheta) {
          const double sintheta = sin(theta);
          const double costheta = cos(theta);
          const double dcostheta = cos(theta - dtheta/2) - cos(theta + dtheta/2);
          /*
            // Integrating around phi is not strictly needed, since
            // the system has a cylindrical symmetry, but could be
            // nice, as it give us an average over grid points.

          const double dphi = 2*M_PI/ceil(delta_r*2*M_PI/dv);
          const double darea = delta_r*delta_r*dcostheta*dphi;
          for (double phi = dphi/2; phi < 2*M_PI; phi += dphi) {
            const Cartesian r1(delta_r*cos(phi)*sintheta,
                               delta_r*sin(phi)*sintheta,
                               z0 + delta_r*costheta);
            double g2 = pairdists[version](gsigma, density, nA, n3, r0, r1);
            da_dz += density(r0)*density(r1)*g2*darea;
          }
          */
          const double darea = delta_r*delta_r*dcostheta*2*M_PI;
          const Cartesian r1(delta_r*sintheta, 0, z0 + delta_r*costheta);
          double g2 = pairdists[version](gsigma, density, nA, n3, r0, r1);
          da_dz += density(r0)*density(r1)*g2*darea;
        }
        fprintf(out, "%g %g\n",z0,da_dz);
      }
      fclose(out);
      char z0_string[50];
      sprintf(z0_string,"%s a1 integral, dv = %g, delta_r = %g",fun[version],dv,delta_r);
      took(z0_string);
    }
  }
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

    FILE *fout = fopen("papers/pair-correlation/figs/wallsfillingfracInfo.txt", "a");
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

  took("Plotting stuff");
}

int main(int, char **) {
  FILE *fout = fopen("papers/pair-correlation/figs/wallsfillingfracInfo.txt", "w");
  fclose(fout);
  read_mc();
  printf("the last g = %g\n", g[10*1300+1300+1]);
  printf("Done with read\n");
  for (double this_eta = 0.1; this_eta < 0.55; this_eta += 0.1) {
    run_walls(this_eta, "WB", WB);
  }
  // Just create this file so make knows we have run.
  if (!fopen("papers/pair-correlation/figs/walls.dat", "w")) {
    printf("Error creating walls.dat!\n");
    return 1;
  }
  fflush(stdout);
  return 0;
}
