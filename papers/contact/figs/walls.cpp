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

// Here we set up the lattice.
static double width = 30;
const double dw = 0.001;
const double spacing = 3; // space on each side

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

Functional WB = HardSpheresNoTensor2(1.0);
Functional WBm2 = WhiteBearMarkII(1.0);
Functional WBT = TensorWhiteBear(1.0);

const int numiters = 25;


void z_plot(const char *fname, const Grid &a, const Grid &b, const Grid &c, const Grid &d,
            const Grid &e, const Grid &f, const Grid &g, const Grid &h) {
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
    double dhere = d(x,y,z);
    double ehere = e(x,y,z);
    double fhere = f(x,y,z);
    double ghere = g(x,y,z);
    double hhere = h(x,y,z);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", here[2],
            ahere, bhere, chere, dhere, ehere, fhere, ghere, hhere);
  }
  fclose(out);
}


void run_walls(double eta, const char *name, Functional fhs) {
  //printf("Filling fraction is %g\n", eta);

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
  //f.run_finite_difference_test("foobar", 1, potential);

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
  
  char *plotname = (char *)malloc(1024);
  Grid gSigmaSGrid(gd, gSigmaS(1.0)(1, gd, density));
  Grid gSigmaAGrid(gd, gSigmaA(1.0)(1, gd, density));
  if (strlen(name) == 4) { 
    gSigmaSGrid = gSigmaSm2(1.0)(1, gd, density);
    gSigmaAGrid = gSigmaAm2(1.0)(1, gd, density);
  }
  Grid gross_correlation(gd, CorrelationGrossCorrect(1.0)(1, gd, density));
  Grid n0(gd, ShellConvolve(1)(1, density)/(4*M_PI));
  Grid nA(gd, ShellConvolve(2)(1, density)/(4*M_PI*4));
  Grid yuwu_correlation(gd, YuWuCorrelation_S(1.0)(1, gd, density));
  sprintf(plotname, "papers/contact/figs/walls%s-%04.2f.dat", name, eta);
  z_plot(plotname, density, density /* fixme hokey */, gSigmaSGrid, yuwu_correlation,
         gSigmaAGrid, n0, gross_correlation, nA);
  free(plotname);

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
    
    FILE *fout = fopen("papers/contact/figs/wallsfillingfracInfo.txt", "a");
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
  FILE *fout = fopen("papers/contact/figs/wallsfillingfracInfo.txt", "w");
  fclose(fout);
  for (double eta = 0.1; eta < 0.6; eta+=0.1) {
    run_walls(eta, "WB", WB);
    run_walls(eta, "WBT", WBT);
    run_walls(eta, "WBm2", WBm2);
  }
  // Just create this file so make knows we have run.
  if (!fopen("papers/contact/figs/walls.dat", "w")) {
    printf("Error creating walls.dat!\n");
    return 1;
  }
  return 0;
}
