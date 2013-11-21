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
#include "new/LJHughesFast.h"
#include <popt.h>

static const double nm = 18.8972613;
static const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin
static const double angstrom = 1.88972613; // An angstrom in atomic units

// Here we set up the lattice.
double zmax = 2.5*nm;
double ymax = 2.5*nm;
double xmax = 2.5*nm;
double dx = 0.1*nm;
int time_me_please = 0;

const char *filename = 0;

const int popt_double = POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT;

int main(int argc, char *argv[]) {
  // ---------------------------------------------------------------------------
  // Set values from parameters
  // ---------------------------------------------------------------------------
  poptOption optionsTable[] = {
    {"xmax", '\0', popt_double, &xmax, 0, "X cell width", "width"},
    {"ymax", '\0', popt_double, &ymax, 0, "Y cell width", "width"},
    {"zmax", '\0', popt_double, &zmax, 0, "Z cell width", "width"},
    {"dx", '\0', popt_double, &dx, 0, "Resolution", "dx"},
    {"time", 't', POPT_ARG_NONE, &time_me_please, 0, "Display timing information", ""},
    {"outfile", 'o', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     "output filename", "filename"},
    POPT_AUTOHELP
    POPT_TABLEEND
  };
  poptContext optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);

  int c = 0;
  // go through arguments, set them based on optionsTable
  while((c = poptGetNextOpt(optCon)) >= 0);
  if (c < -1) {
    fprintf(stderr, "\n%s: %s\n", poptBadOption(optCon, 0), poptStrerror(c));
    return 1;
  }
  poptFreeContext(optCon);

  LJHughes f(xmax, ymax, zmax, dx);
  f.R() = 3.03420*angstrom/2; // hard sphere radius of water in bohr from Clark et al, table 6
  f.kT() = 298.15*kB; // room temperature in Hartree

  f.epsilon_association() = 1400.00*kB; // epsilonAB from Clark et al, table 6
  f.kappa_association() = 0.0381876*(3.03420*angstrom*3.03420*angstrom*3.03420*angstrom); // kappaAB from Clark et al code

  f.epsilon_dispersion() = 250.0*kB; // epsilon_dispersion from Clark et al, table 6
  f.lambda_dispersion() = 1.78890; // lambda_dispersion from Clark et al, table 6
  f.length_scaling() = 0.353; // length-scaling parameter that DJR added to fix surface tension

  f.ljepsilon() = 0.1;
  f.ljsigma() = 0.1;

  f.mu() = 0;

  // Initialize the potential to zero...
  f.Veff() = 0.0;

  Vector n = f.get_n();
  Vector rx = f.get_rx();
  Vector ry = f.get_ry();
  Vector rz = f.get_rz();
  FILE *o = fopen(filename, "w");
  if (!o) {
    printf("Error opening file %s\n", filename);
    exit(1);
  }
  for (int i=0; i<f.Nx(); i++) {
    fprintf(o, "%g\t%g\t%g\t%g\n", rx[i], ry[i], rz[i], n[i]);
  }
  fclose(o);
}
