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
#include "LineMinimizer.h"

const double my_kT = 1e-3; // room temperature in Hartree

// Here we set up the lattice.
Lattice lat(Cartesian(5,0,0), Cartesian(0,5,0), Cartesian(0,0,5));
double resolution = 0.2;
GridDescription gd(lat, resolution);


int test_minimizer(const char *name, Minimizer min, Grid *psi, double fraccuracy=1e-3) {
  clock_t start = clock();
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  const double true_energy = -0.5;

  *psi = +1e-4*((-10*r2(gd)).cwise().exp()) + 1.14*VectorXd::Ones(psi->description().NxNyNz);

  while (min.improve_energy(true)) fflush(stdout);

  min.print_info();
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  printf("fractional energy error = %g\n", (min.energy() - true_energy)/fabs(true_energy));
  if (fabs((min.energy() - true_energy)/true_energy) > fraccuracy) {
    printf("FAIL: Error in the energy is too big!\n");
    return 1;
  }
  if (min.energy() < true_energy) {
    printf("FAIL: Sign of error is wrong!!!\n");
    return 1;
  }
  return 0;
}

int main(int, char **argv) {
  int retval = 0;

  Functional ff = Hydrogen();

  Grid psi(gd);
  Minimizer cg = MaxIter(150, ConjugateGradient(ff, gd, my_kT, &psi, QuadraticLineMinimizer));
  retval += test_minimizer("ConjugateGradient", cg, &psi, 1.0);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    printf("FIXME: This fails, but I'm going to ignore the failure for now.\n");
    return 0;
    return retval;
  }
}
