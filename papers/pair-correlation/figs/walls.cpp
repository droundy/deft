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
#include "errno.h"
#include "sys/stat.h" // for mkdir

double radial_distribution(double gsigma, double r);
double py_rdf (double eta, double r);

// Maximum and spacing values for plotting
const double zmax = 20;
const double xmax = 10;
const double dx = 0.1;


// The functions for different ways of computing the pair distribution function.
double pairdist_simple(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) {
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  return (radial_distribution(gsigma(r0), r) + radial_distribution(gsigma(r1), r))/2;
}
double pairdist_nA(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) {
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  return (radial_distribution(gsigma(r0), r)/nA(r0) + radial_distribution(gsigma(r1), r)/nA(r1))
    /(1/nA(r0) + 1/nA(r1));
}
double pairdist_py(const Grid &gsigma, const Grid &n, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) {
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  const double eta = 4/3*M_PI*1*1*1*(n(r0) + n(r1))/2;
  return py_rdf(eta, r);
}
double pairdist_fischer(const Grid &gsigma, const Grid &n, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) {
  // This implements the pair distribution function of Fischer and
  // Methfessel from the 1980 paper.  The py_rdf below should be the
  // true radial distribution function for a homogeneous hard-sphere
  // fluid with packing fraction eta.
  const Cartesian r01 = Cartesian(r0 - r1);
  const double r = sqrt(r01.dot(r01));
  const double eta = n3(Cartesian(0.5*(r0+r1)));
  return py_rdf(eta, r);
}

const char *fun[] = {
  "simple",
  "nA",
  "py",
  "fischer"
};
double (*pairdists[])(const Grid &gsigma, const Grid &density, const Grid &nA, const Grid &n3, Cartesian r0, Cartesian r1) = {
  pairdist_simple,
  pairdist_nA,
  pairdist_py,
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
  double d = r/2 - 1;
  if (d < 0)
    return 0;
  return 1 + (gsigma-1)*exp(-a0*d) + a1*(gsigma-1)*sin(a2*d)*exp(-a3*d)
           - a4*(gsigma-1)*(gsigma-1)*sin(a5*d)*exp(-a6*d);
}

double py_rdf (double eta, double r) {
  // taken from the paper by Trokhymchuk et. al.
  const double sigma = 1.0;
  const double gsigma = 1.0/4.0/eta*((1.0 + eta + eta*eta - 2.0/3.0*eta*eta*eta
                      - 2.0/3.0*eta*eta*eta*eta)/(1.0 - eta)/(1.0 - eta)/(1.0 - eta) - 1.0);
  const double rstar = sigma*(2.0116 - 1.0647*eta + 0.0538*eta*eta);
  const double d = pow((2.0*eta*(eta*eta - 3.0*eta - 3.0 + sqrt(3.0*(eta*eta*eta*eta
	               - 2.0*eta*eta*eta + eta*eta + 6.0*eta + 3.0)))), 1.0/3.0);
  const double g_m = 1.0286 - 0.6095*eta + 3.5781*eta*eta - 21.3651*eta*eta*eta
                   + 42.6344*eta*eta*eta*eta - 33.8485*eta*eta*eta*eta*eta;
  const double alpha = (44.554 + 79.868*eta + 116.432*eta*eta - 44.652*exp(2.0*eta))/sigma;
  const double alpha0 = 2.0*eta/(1.0 - eta)*(-1.0 + d/4.0/eta - eta/2.0/d)/sigma;
  const double beta = (-5.022 + 5.857*eta + 5.089*exp(-4.0*eta))/sigma;
  const double beta0 = 2.0*eta/(1.0-eta)*sqrt(3.0)*(-d/4.0/eta - eta/2.0/d)/sigma;
  const double mu = (2.0*eta/(1.0 - eta)*(-1.0 - d/2.0/eta - eta/d))/sigma;
  const double gamma = atan(-sigma/beta0*((alpha0*sigma*(alpha0*alpha0 + beta0*beta0)
                     - mu*sigma*(alpha0*alpha0 + beta0*beta0))*(1.0 + eta/2.0)
                     + (alpha0*alpha0 + beta0*beta0 - mu*alpha0)*(1.0 + 2.0*eta)));
  const double kappa = (4.674*exp(-3.935*eta) + 3.536*exp(-56.270*eta))/sigma;
  const double omega = (-0.682*exp(-24.697*eta) + 4.720 + 4.450*eta)/sigma;
  const double delta = -omega*rstar - atan((kappa*rstar + 1.0)/omega/rstar);
  const double C = rstar*(g_m - 1.0)*exp(kappa*rstar)/cos(omega*rstar + delta);
  const double B = (g_m - (sigma*gsigma/rstar)*exp(mu*(rstar - sigma)))
                 / (cos(beta*(rstar - sigma) + gamma)*exp(alpha*(rstar-sigma))
		    - cos(gamma)*exp(mu*(rstar - sigma)))*rstar;
  const double A = sigma*gsigma - B*cos(gamma);
  if (r < sigma)
    return 0;
  else if ( r < rstar) {
   const double g_dep = A/r*exp(mu*(r-sigma))
                      + B/r*cos(beta*(r-sigma) + gamma)*exp(alpha*(r-sigma));
   return g_dep;
  }
  else {
    const double g_str = 1 + C/r*cos(omega*r + delta)*exp(-kappa*r);
    return g_str;
  }
}

double notinwall(Cartesian r) {
  const double z = r.z();
  if (fabs(z) > spacing) {
      return 1;
  }
  return 0;
}

static void took(const char *name) {
  assert(name); // so it'll count as being used...
  // static clock_t last_time = clock();
  // clock_t t = clock();
  // double peak = peak_memory()/1024.0/1024;
  // printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  // fflush(stdout);
  // last_time = t;
}

Functional WB = HardSpheresNoTensor(1.0);
//Functional WBm2 = HardSpheresWBm2(1.0);
//Functional WBT = HardSpheresWBFast(1.0);

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
  // hypothetical loop
  //for (double z=0; z<gd.Lat.a3().z(); z+= 0.015) {
  //  double ahere = a(Cartesian(0,0,z));
  //  printf("a(%g) = %g\n", z, ahere);
  //}
  fclose(out);
}


void run_walls(double eta, const char *name, Functional fhs) {
  // Generates a data file for the pair distribution function, for filling fraction eta
  // and distance of first sphere from wall of z0. Data saved in a table such that the
  // columns are x values and rows are z1 values.

  Functional f = OfEffectivePotential(fhs + IdealGas());
  double mu = find_chemical_potential(f, 1, eta/(4*M_PI/3));
  f = OfEffectivePotential(fhs + IdealGas()
                           + ChemicalPotential(mu));

  const double zmax = width + 2*spacing;
  Lattice lat(Cartesian(dw,0,0), Cartesian(0,dw,0), Cartesian(0,0,zmax));
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
  Grid gsigma(gd, Correlation_A2(1.0)(1, gd, density));
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
  for (double z0 = 3.05; z0 < 13; z0 += .1) {
    // For each z0, we now pick one of our methods for computing the
    // pair distribution function:
    for (int version = 0; version < numplots; version++) {
      sprintf(plotname, "papers/pair-correlation/figs/walls/walls%s-%s-pair-%04.2f-%1.2f.dat",
              name, fun[version], eta, z0);
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

  for (double eta = 0.1; eta < 0.6; eta += 0.1) {
    run_walls(eta, "WB", WB);
    //run_walls(eta, "WBT", WBT);
    //run_walls(eta, "WBm2", WBm2);
  }
  // Just create this file so make knows we have run.
  if (!fopen("papers/pair-correlation/figs/walls.dat", "w")) {
    printf("Error creating walls.dat!\n");
    return 1;
  }
  return 0;
}
