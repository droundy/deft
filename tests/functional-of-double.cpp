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
#include "Functionals.h"

int test_functional(const char *name, Functional f, double n, double fraccuracy=1e-14) {
  printf("\n**************************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s of %10g *\n", name, n);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("**************************\n\n");

  // Here we set up the lattice.
  Lattice lat(Cartesian(4.2,0,-0.1), Cartesian(0.1,5.2,0), Cartesian(0,0,5));
  double resolution = 0.2;
  GridDescription gd(lat, resolution);
  Grid nr(gd, n*VectorXd::Ones(gd.NxNyNz));

  const double Edouble = f(n)*gd.Lat.volume();
  const double Egrid = f(nr);

  printf("fractional error = %g\n", (Edouble - Egrid)/fabs(Edouble));
  if (fabs((Edouble - Egrid)/Edouble) > fraccuracy) {
    printf("Error in the energy is too big!\n");
    return 1;
  }
  return 0;
}

int main(int, char **argv) {
  int retval = 0;
  const double kT = 1e-3;
  const FieldFunctional n = EffectivePotentialToDensity(kT);

  {
    Functional attr = GaussianPolynomial(-0.32, 0.5, 2);
    retval += test_functional("Attractive Gaussian", attr, 0.1, 1e-13);
    Functional repul = GaussianPolynomial(0.32, 0.25, 4);
    retval += test_functional("Repulsive Gaussian", repul, 0.1, 1e-12);
    retval += test_functional("sum of gaussians", attr + repul, 0.1, 1e-13);
    retval += test_functional("other sum of gaussians", repul + attr, 0.1, 1e-13);
  }

  {
    Functional f = IdealGas(kT);
    retval += test_functional("Ideal gas", f, 1e-9, 1e-13);
    retval += test_functional("Ideal gas", f, 1e-3, 1e-12);
    retval += test_functional("Ideal gas of V", f(n), -kT*log(1e-9), 1e-13);
    retval += test_functional("Ideal gas of V", f(n), -kT*log(1e-3), 1e-12);
  }

  {
    Functional f = ChemicalPotential(0.1);
    retval += test_functional("chemical potential", f, 1e-9, 1e-12);
    retval += test_functional("chemical potential", f, 1e9, 1e-14);
    retval += test_functional("chemical potential", f, 1e-2, 1e-13);
    retval += test_functional("chemical potential of V", f(n), -kT*log(1e-9), 1e-12);
    retval += test_functional("chemical potential of V", f(n), -kT*log(1e-3), 1e-12);
  }

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
