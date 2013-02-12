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

// Maximum and spacing values for plotting, saved for use by plot-walls.py
const double zmax = 20;
const double xmax = 10;
const double dx = 0.1;

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
  if (d <= 0)
    return 0;
  return 1 + (gsigma-1)*exp(-a0*d) + a1*(gsigma-1)*sin(a2*d)*exp(-a3*d)
           - a4*(gsigma-1)*(gsigma-1)*sin(a5*d)*exp(-a6*d);
}

double py_rdf (double eta, double r) {
  // taken from the paper by Trokhymchuk et. al.
  const double sigma = 1;
  const double gsigma = 1/4/eta*((1 + eta + eta*eta - 2/3*eta*eta*eta
                      - 2/3*eta*eta*eta*eta)/(1 - eta)/(1 - eta)/(1 - eta) - 1);
  const double rstar = sigma*(2.0116 - 1.0647*eta + 0.0538*eta*eta);
  const double d = pow((2*eta*(eta*eta - 3*eta - 3 + sqrt(3*(eta*eta*eta*eta
	               - 2*eta*eta*eta + eta*eta + 6*eta + 3)))), 1.0/3.0);
  const double g_m = 1.0286 - 0.6095*eta + 3.5781*eta*eta - 21.3651*eta*eta*eta
                   + 42.6344*eta*eta*eta*eta - 33.8485*eta*eta*eta*eta*eta;
  const double alpha = (44.554 + 79.868*eta + 116.432*eta*eta - 44.652*exp(2*eta))/sigma;
  const double alpha0 = 2*eta/(1 - eta)*(-1 + d/4/eta - eta/2/d)/sigma;
  const double beta = (-5.022 + 5.857*eta + 5.089*exp(-4*eta))/sigma;
  const double beta0 = 2*eta/(1-eta)*sqrt(3)*(-d/4/eta - eta/2/d)/sigma;
  const double mu = (2*eta/(1 - eta)*(-1 - d/2/eta - eta/d))/sigma;
  const double gamma = atan(-sigma/beta0*((alpha0*sigma*(alpha0*alpha0 + beta0*beta0)
                     - mu*sigma*(alpha0*alpha0 + beta0*beta0))*(1 + eta/2)
                     + (alpha0*alpha0 + beta0*beta0 - mu*alpha0)*(1 + 2*eta)));
  const double kappa = (4.674*exp(-3.935*eta) + 3.536*exp(-56.270*eta))/sigma;
  const double omega = (-0.682*exp(-24.697*eta) + 4.720 + 4.450*eta)/sigma;
  const double delta = -omega*rstar - atan((kappa*rstar + 1)/omega/rstar);
  const double C = rstar*(g_m - 1)*exp(kappa*rstar)/cos(omega*rstar + delta);
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
  //static clock_t last_time = clock();
  //clock_t t = clock();
  assert(name); // so it'll count as being used...
  //double peak = peak_memory()/1024.0/1024;
  //printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  //last_time = t;
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

void plot_pair_distribution(const char *fname, const char *fname2, double z0,
                            const Grid &gsigma, const Grid &density, const Grid &nA) {
  FILE *out = fopen(fname, "w");
  FILE *out2 = fopen(fname2, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    return;
  }
  if (!out2) {
    fprintf(stderr, "Unable to create file %s!\n", fname2);
    return;
  }

  // the +1 for z0 and z1 are to shift the plot over, so that a sphere touching the wall
  // is at z = 0, to match with the monte carlo data
  z0 += 1;
  double gsigma0 = gsigma(Cartesian(0, 0, z0));
  double density0 = density(Cartesian(0, 0, z0+2)); //density starts at z=3 for some reason
  double nA0 = nA(Cartesian(0, 0, z0));
  const GridDescription gd = gsigma.description();
  for (double x = 0; x < xmax - dx/2; x += dx) {
    for (double z1 = 1; z1 < zmax + 1 - dx/2; z1 += dx) {
      double gsigma1 = gsigma(Cartesian(0, 0, z1));
      double density1 = density(Cartesian(0, 0, z1+2));
      double nA1 = nA(Cartesian(0, 0, z1));
      double r = sqrt((z1 - z0)*(z1 - z0) + x*x);
      double g2 = 0;
      double g2b = 0;
      if (r >= 2) {
        //g2 = (radial_distribution(gsigma0, r) + radial_distribution(gsigma1, r))/2;
        //g2b = (radial_distribution(gsigma0, r)*density2 + radial_distribution(gsigma2, z1+2)*density0)/(density0+density2);
        //g2 = density1;//(radial_distribution(gsigma0, r));//*density1)/(density0 + density1);
        //g2b = (radial_distribution(gsigma1, r));//*density0)/(density0 + density1);
        //g2 = radial_distribution(gsigma0, r) + radial_distribution(gsigma1, z1+1);
        //g2b = (nA1*radial_distribution(gsigma0, r) + nA0*radial_distribution(gsigma1, r))
	/// (nA0 + nA1);
        //        g2 = (radial_distribution(gsigma0, r)*density1 + radial_distribution(gsigma1, r)*density0)/(density0 + density1);
        // g2 = radial_distribution(gsigma0, r);
        //   g2b = radial_distribution(gsigma1, r)*density0/density1;

        //g2 = radial_distribution(gsigma0, r)*nA1/(nA0 + nA1) + radial_distribution(gsigma1, r)*nA0/(nA0 + nA1);
        //g2b = (radial_distribution(gsigma1, r)/density1 + radial_distribution(gsigma0, r)/density0)*(density0);
        //g2 = radial_distribution(gsigma1, r)/density1*(density0 + density1) + radial_distribution(gsigma0, r)/density0*(density0 + density1);
        g2 = (radial_distribution(gsigma0, r)/nA1 + radial_distribution(gsigma1, r)/nA0)
          /(1/nA0 + 1/nA1);
        g2b = (radial_distribution(gsigma0, r)/density0
             + radial_distribution(gsigma1, r)/density1)
          /(1/density0 + 1/density1);
        const double a = radial_distribution(gsigma1, r)/density1/(1/density0 + 1/density1);
        const double b = radial_distribution(gsigma0, r)/density0/(1/density0 + 1/density1);
        g2 = a + b;
      }
      fprintf(out, "%g\t", g2);
      fprintf(out2, "%g\t", g2b);
    }
    fprintf(out, "\n");
    fprintf(out2, "\n");
  }
  fclose(out);
  fclose(out2);
}

void plot_py_rdf(const char *fname, double z0, const Grid &gsigma, double eta) {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unabme to create file %s!\n", fname);
    return;
  }
  double gsigma0 = gsigma(Cartesian(0, 0, z0+1));
  const GridDescription gd = gsigma.description();
  for (double x = 0; x < xmax - dx/2; x += dx) {
    for (double z1 = 0; z1 < zmax - dx/2; z1 += dx) {
      double gsigma1 = gsigma(Cartesian(0, 0, z1+1));
      double r = sqrt((z0-z1)*(z0-z1) + x*x);
      double g2py = 0;
      if (r > 2) {
	g2py = py_rdf(eta, r/2);
      }
      fprintf(out, "%g\t", g2py);
    }
    fprintf(out, "\n");
  }
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
  char *plotname2 = new char[1024];
  char *plotname3 = new char[1024];
  Grid gsigma(gd, Correlation_A2(1.0)(1, gd, density));
  Grid nA(gd, ShellConvolve(2)(1, density)/(4*M_PI*4));

  sprintf(plotname, "papers/pair-correlation/figs/walls%s-%04.2f.dat", name, eta);
  z_plot(plotname, density, gsigma, nA);

  // here you choose the values of z0 to use
  for (double z0 = 0; z0 < 3; z0 += .5) {
    sprintf(plotname, "papers/pair-correlation/figs/walls%s-pair-%04.2f-%g.dat",
	    name, eta, z0);
    sprintf(plotname2, "papers/pair-correlation/figs/walls%s-pair-%04.2f-%g-b.dat",
	    name, eta, z0);
    sprintf(plotname3, "papers/pair-correlation/figs/walls%s-pair-py-%04.2f-%g.dat",
	    name, eta, z0);
    plot_pair_distribution(plotname, plotname2, z0, gsigma, density, nA);
    plot_py_rdf(plotname3, z0, gsigma, eta);
  }
  delete[] plotname;
  delete[] plotname2;
  delete[] plotname3;
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
  FILE *out = fopen("papers/pair-correlation/figs/constants.dat", "w");
  if (!out) {
    fprintf(stderr, "Unable to create file constants.dat!\n");
    return 1;
  }
  fprintf(out, "%g\t%g\t%g", zmax, xmax, dx);
  fclose(out);

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
