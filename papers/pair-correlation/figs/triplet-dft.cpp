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

//const double dx = 0.1;
const double dx = .1/1.5;//This is a bit smaller than 1/sqrt(2),
                          //which allows for a speration delta of 1
                          //when going around the circular wall and
                          //not having to worry about a gridpoint
                          //being within the wall.

double mc_r_step;
const int num_eta = 10;
const double eta_step = 0.5/num_eta;
const int num_r_in_gmc = 1300;
const int gsize = num_r_in_gmc*(num_eta+1);
double * g = new double[gsize];

// The functions for different ways of computing the pair distribution function.

double gsigma_to_eta(const double gs) {
  if (gs <= 1) return 0;
  const double c = -pow(fabs(sqrt(3.0)*sqrt(27.0*gs*gs*gs*gs - 2.0*gs*gs*gs)-9.0*gs*gs), 1.0/3.0);
  return c/(pow(6.0, 2.0/3.0)*gs) + 1.0/(pow(6.0, 1.0/3.0)*c) + 1.0;
}

double pairdist_this_work(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3,
                        const Grid &nbar_sokolowski, Cartesian r0, Cartesian r1) {
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  return (radial_distribution(gsigma(r0), r) + radial_distribution(gsigma(r1), r))/2;
}

double pairdist_this_work_mc(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3,
                             const Grid &nbar_sokolowski,  Cartesian r0, Cartesian r1) {
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  const double eta0 = gsigma_to_eta(gsigma(r0));
  const double eta1 = gsigma_to_eta(gsigma(r1));
  return (mc(eta0, r, mc_r_step, g) + mc(eta1, r, mc_r_step, g))/2;
}
double pairdist_fischer(const Grid &gsigma, const Grid &n, const Grid &nA, const Grid &n3,
                        const Grid &nbar_sokolowski, Cartesian r0, Cartesian r1) {
  // This implements the pair distribution function of Fischer and
  // Methfessel from the 1980 paper.  The py_rdf below should be the
  // true radial distribution function for a homogeneous hard-sphere
  // fluid with packing fraction eta.
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  const double eta = n3(Cartesian(0.5*(r0+r1)));
  return mc(eta, r, mc_r_step, g);
}
double pairdist_sokolowski(const Grid &gsigma, const Grid &n, const Grid &nA, const Grid &n3,
                           const Grid &nbar_sokolowski, Cartesian r0, Cartesian r1) {
  // This implements the pair distribution function of Sokolowski and
  // Fischer from the 1992 paper.
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));

  const double eta0 = nbar_sokolowski(r0)*(4.0/3.0*M_PI);
  const double eta1 = nbar_sokolowski(r1)*(4.0/3.0*M_PI);
  const double eta = (eta0 + eta1)/2.0;

  return mc(eta, r, mc_r_step, g);
}

const char *fun[] = {
  "this-work",
  "this-work-mc",
  "sokolowski",
  "fischer"
};
double (*pairdists[])(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3,
                      const Grid &nbar_sokolowski, Cartesian r0, Cartesian r1) = {
  pairdist_this_work,
  pairdist_this_work_mc,
  pairdist_sokolowski,
  pairdist_fischer
};
const int numplots = sizeof fun/sizeof fun[0];


// Here we set up the lattice.
static double width = 30;
const double dw = 0.001;

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
  char *fname = new char[4096];
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
  for (int i=0; i<num_r_in_gmc; i++) {
    g[i] = 1.0;
  }
  for (double eta = eta_step; eta < 0.5001; eta += eta_step) {
    char *fname = new char[4096];
    sprintf(fname, "papers/pair-correlation/figs/gr-%2.2f.dat", eta);
    FILE *in = fopen(fname, "r");
    if (!in) {
      fprintf(stderr, "Unable to open file %s!!!!\n", fname);
      delete[] fname;
      exit(1);
    }
    int i=0;
    while (fscanf(in, " %*g %lg", &num1) == 1) {
      g[int(floor(num_r_in_gmc*(eta/eta_step)))+i] = num1/eta;
      i++;
    }
    if (i!=num_r_in_gmc) {
      printf("There  are not num_r_in_gmc lines in papers/pair-correlation/figs/gr-%2.2f.dat which is a problem in walls.cpp\n", eta);
      exit(1);
    }
    delete[] fname;
    fclose(in);
  }
}

double mc(double local_eta, double r, double r_step, double g[]) {
  int low_eta_i = floor(local_eta/eta_step);
  int high_eta_i = low_eta_i + 1;
  double g_low_eta;
  double g_high_eta;
  double fac;
  int j=0;
  double r_floor;
  if (r < 2) {
    return 0;
  }
  if (high_eta_i > num_eta) {
    //fprintf(stderr, "Error. Trying to interpolate to eta = %g, but our max is %g\n", local_eta, num_eta*eta_step);
    return NAN;
  }
  if (r < 2 +(r_step/2.0)) {
    fac = (r-2.0)/(0.5*r_step);
    g_low_eta = (1-fac)*g[low_eta_i*num_r_in_gmc]+fac*g[low_eta_i*num_r_in_gmc + 1];
    g_high_eta = (1-fac)*g[high_eta_i*num_r_in_gmc]+fac*g[high_eta_i*num_r_in_gmc + 1];
  } else if (r > 14.985) {
    return 1.0;
  } else {
    j = floor((r-2.0+0.5*r_step)/r_step);
    r_floor = 2+(j-0.5)*r_step;
    fac = (r-r_floor)/r_step;
    g_low_eta = (1-fac)*g[low_eta_i*num_r_in_gmc+j] + fac*g[low_eta_i*num_r_in_gmc+j+1];
    g_high_eta = (1-fac)*g[high_eta_i*num_r_in_gmc+j] + fac*g[high_eta_i*num_r_in_gmc+j+1];
  }
  //printf("low_eta_i = %d and high_eta_i = %d\n", low_eta_i*num_r_in_gmc+j, high_eta_i*num_r_in_gmc+j);
  //fflush(stdout);
  fac = (local_eta-(low_eta_i)*eta_step)/eta_step;
  double ghere = (1-fac)*g_low_eta + fac*g_high_eta;
  return (1-fac)*g_low_eta + fac*g_high_eta;
}

double notinsphere(Cartesian r) {
  const double x = r.x();
  const double y = r.y();
  const double z = r.z();
  const double dis = sqrt(z*z+y*y+x*x);
  if (dis<2) {
    return 0;
  }
  return 1;
}

Functional WB = HardSpheresNoTensor2(1.0);

static void took(const char *name) {
  // assert(name); // so it'll count as being used...
  // static clock_t last_time = clock();
  // clock_t t = clock();
  // double peak = peak_memory()/1024.0/1024;
  // double seconds = (t-last_time)/double(CLOCKS_PER_SEC);
  // if (seconds > 120) {
  //   printf("\t\t%s took %.0f minutes and %g M memory\n", name, seconds/60, peak);
  // } else {
  //   printf("\t\t%s took %g seconds and %g M memory\n", name, seconds, peak);
  // }
  // fflush(stdout);
  // last_time = t;
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
  Lattice lat(Cartesian(width,0,0), Cartesian(0,width,0), Cartesian(0,0,width));
  GridDescription gd(lat, 0.1);
  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinsphere);
  f = constrain(constraint, f);
  potential = (eta*constraint + 1e-4*eta*VectorXd::Ones(gd.NxNyNz))/(4*M_PI/3);
  potential = -potential.cwise().log();

  const double approx_energy = (fhs + IdealGas() + ChemicalPotential(mu))(1, eta/(4*M_PI/3))*dw*dw*width;
  const double precision = fabs(approx_energy*1e-1);//1e-4);
  //printf("Minimizing to %g absolute precision...\n", precision);
  { // Put mimizer in block so as to free it when we finish minimizing to save memory.
    Minimizer min = Precision(precision,
                              PreconditionedConjugateGradient(f, gd, 1,
                                                              &potential,
                                                              QuadraticLineMinimizer));
    for (int i=0;min.improve_energy(true) && i<100;i++) {
      double peak = peak_memory()/1024.0/1024;
      double current = current_memory()/1024.0/1024;
      printf("Peak memory use is %g M (current is %g M)\n", peak, current);
      fflush(stdout);
    }
    took("Doing the minimization");
  }
  Grid density(gd, EffectivePotentialToDensity()(1, gd, potential));
  //printf("# per area is %g at filling fraction %g\n", density.sum()*gd.dvolume/dw/dw, eta);
  Grid gsigma(gd, gSigmaA(1.0)(1, gd, density));
  Grid nA(gd, ShellConvolve(2)(1, density)/(4*M_PI*4));
  Grid n3(gd, StepConvolve(1)(1, density));
  Grid nbar_sokolowski(gd, StepConvolve(1.6)(1, density));
  nbar_sokolowski /= (4.0/3.0*M_PI*ipow(1.6, 3));
  // Create the walls directory if it doesn't exist.
  if (mkdir("papers/pair-correlation/figs/walls", 0777) != 0 && errno != EEXIST) {
    // We failed to create the directory, and it doesn't exist.
    printf("Failed to create papers/pair-correlation/figs/walls: %s",
           strerror(errno));
    exit(1); // fail immediately with error code
  }

  // here you choose the values of z0 to use
  //dx is set at beggining of file
  char *plotname = new char[4096];
  for (double z0 = 0; z0 < width/2.0; z0 += dx) {
    // For each z0, we now pick one of our methods for computing the
    // pair distribution function:
    for (int version = 0; version < numplots; version++) {
      sprintf(plotname,
              "papers/pair-correlation/figs/triplet%s-%s-%04.2f-%1.2f.dat",
              name, fun[version], eta, z0);
      FILE *out = fopen(plotname,"w");
      // the +1 for z0 and z1 are to shift the plot over, so that a sphere touching the wall
      // is at z = 0, to match with the monte carlo data
      const Cartesian r0(0,0,z0);
      for (double x = 0; x < width/2; x += dx) {
        for (double z1 = -width/2.0; z1 < width/2.0; z1 += dx) {
          const Cartesian r1(x,0,z1);
          double g2 = pairdists[version](gsigma, density, nA, n3, nbar_sokolowski, r0, r1);
          double n_bulk = (3.0/4.0/M_PI)*eta;
          double g3 = g2*density(r0)*density(r1)/n_bulk/n_bulk;
          fprintf(out, "%g\t", g3);
        }
        fprintf(out, "\n");
      }
      fclose(out);
    }
  }
  delete[] plotname;
  took("Dumping the triplet dist plots");
  char *plotname_path = new char[4096];
  for (int version = 0; version < numplots; version++) {
    sprintf(plotname_path,
            "papers/pair-correlation/figs/triplet%s-path-%s-%04.2f.dat",
            name, fun[version], eta);
    FILE *out_path = fopen(plotname_path, "w");
    if (!out_path) {
      fprintf(stderr, "Unable to create file %s!\n", plotname_path);
      return;
    }
    fprintf(out_path, "# unused\tg3\tz\tx\n");

    const double delta = .1; //this is the value of radius of the
                              //particle as it moves around the
                              //contact sphere on its path
    int num = 100; //This is the same num that is in plot-path.py,
                  //splits up the theta part of path just like there
    const Cartesian r0(0,0, 2.0+delta);
    const double max_theta = M_PI*2.0/3;
    const double ds = 0.001; // step size to use in path plots
    for (double z = 7; z >= 3*(2.0 + delta); z-=ds) {
      const Cartesian r1(0,0,z);
      double g2_path = pairdists[version](gsigma, density, nA, n3, nbar_sokolowski, r0, r1);
      double n_bulk = (3.0/4.0/M_PI)*eta;
      double g3 = g2_path*density(r0)*density(r1)/n_bulk/n_bulk;
      fprintf(out_path,"0\t%g\t%g\t%g\n", g3, r1[2], r1[0]);
    }
    const double dtheta = ds/2;
    for (double theta = 0; theta <= max_theta; theta += dtheta){
      const Cartesian r1((2.0+delta)*sin(theta), 0, (2.0+delta)*(2+cos(theta)));
      double g2_path = pairdists[version](gsigma, density, nA, n3, nbar_sokolowski, r0, r1);
      double n_bulk = (3.0/4.0/M_PI)*eta;
      double g3 = g2_path*density(r0)*density(r1)/n_bulk/n_bulk;
      fprintf(out_path,"0\t%g\t%g\t%g\n", g3, r1[2], r1[0]);
    }
    for (double x = (2.0+delta)*sqrt(3)/2; x<=6; x+=ds){
      const Cartesian r1(x, 0, 1.0+delta/2);
      double g2_path = pairdists[version](gsigma, density, nA, n3, nbar_sokolowski, r0, r1);
      double n_bulk = (3.0/4.0/M_PI)*eta;
      double g3 = g2_path*density(r0)*density(r1)/n_bulk/n_bulk;
      fprintf(out_path,"0\t%g\t%g\t%g\n", g3, r1[2], r1[0]);
    }
    fclose(out_path);
  }
  for (int version = 0; version < numplots; version++) {
    sprintf(plotname_path,
            "papers/pair-correlation/figs/triplet-path-inbetween-%s-%04.2f.dat",
            fun[version], eta);
    FILE *out_path = fopen(plotname_path, "w");
    if (!out_path) {
      fprintf(stderr, "Unable to create file %s!\n", plotname_path);
      return;
    }
    fprintf(out_path, "# unused\tg3\tz\tx\n");

    const double delta = .1; //this is the value of radius of the
                              //particle as it moves around the
                              //contact sphere on its path
    int num = 100; //This is the same num that is in plot-path.py,
                  //splits up the theta part of path just like there
    const Cartesian r0(0,0, 4.0+2*delta);
    const double max_theta = M_PI;
    const double ds = 0.001; // step size to use in path plots
    for (double z = 10; z >= 3*(2.0 + delta); z-=ds) {
      const Cartesian r1(0,0,z);
      double g2_path = pairdists[version](gsigma, density, nA, n3, nbar_sokolowski, r0, r1);
      double n_bulk = (3.0/4.0/M_PI)*eta;
      double g3 = g2_path*density(r0)*density(r1)/n_bulk/n_bulk;
      fprintf(out_path,"0\t%g\t%g\t%g\n", g3, r1[2], r1[0]);
    }
    const double dtheta = ds/2;
    for (double theta = 0; theta <= max_theta; theta += dtheta){
      const Cartesian r1((2.0+delta)*sin(theta), 0, (2.0+delta)*(1+cos(theta)));
      double g2_path = pairdists[version](gsigma, density, nA, n3, nbar_sokolowski, r0, r1);
      double n_bulk = (3.0/4.0/M_PI)*eta;
      double g3 = g2_path*density(r0)*density(r1)/n_bulk/n_bulk;
      fprintf(out_path,"0\t%g\t%g\t%g\n", g3, r1[2], r1[0]);
    }
    for (double x = 0; x<=6; x+=ds){
      const Cartesian r1(x, 0, 2.0+delta);
      double g2_path = pairdists[version](gsigma, density, nA, n3, nbar_sokolowski, r0, r1);
      double n_bulk = (3.0/4.0/M_PI)*eta;
      double g3 = g2_path*density(r0)*density(r1)/n_bulk/n_bulk;
      fprintf(out_path,"0\t%g\t%g\t%g\n", g3, r1[2], r1[0]);
    }
    fclose(out_path);
  }
  delete[] plotname_path;
}

int main(int, char **) {
  read_mc();
  for (double this_eta = 0.3; this_eta < 0.45; this_eta += 0.1) {
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

