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
#include "OptimizedFunctionals.h"
#include "equation-of-state.h"

int retval = 0;

Functional SoftFluid(double sigma, double epsilon, double mu);

// Here we set up the lattice.
static double width = 15;
const double dx = 0.001;
const double dw = 0.001;
const double spacing = 1.5; // space on each side

double notinwall(Cartesian r) {
  const double z = r.z();
  if (fabs(z) > spacing) {
      return 1;
  }
  return 0;
}

int run_walls(double eta, double temp) {
  Functional f = OfEffectivePotential(SoftFluid(1, 1, 0));
  const double mu = find_chemical_potential(f, temp, eta/(4*M_PI/3));
  f = OfEffectivePotential(SoftFluid(1, 1, mu));

  const double zmax = width + 2*spacing;
  Lattice lat(Cartesian(dw,0,0), Cartesian(0,dw,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, dx);

  Grid constraint(gd);
  constraint.Set(notinwall);
  //f = constrain(constraint, f);

  Grid potential(gd);
  potential = (eta*constraint + 1e-4*eta*VectorXd::Ones(gd.NxNyNz))/(4*M_PI/3);
  potential = -temp*potential.cwise().log();

  char buffer[1024];
  sprintf(buffer, "eta=%g   kT/V0=%g", eta, temp);

  return f.run_finite_difference_test(buffer, temp, potential, NULL);
}

int main(int, char **argv) {
  int retval = 0;

  retval += run_walls(0.4, 0.01);
  //retval += run_walls(0.3, 0.01);
  //retval += run_walls(0.2, 0.01);
  //retval += run_walls(0.1, 0.01);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
