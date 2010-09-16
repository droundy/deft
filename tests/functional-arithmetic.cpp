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

int test_functionals(const char *name, FieldFunctional f1, FieldFunctional f2, double n, double fraccuracy=1e-14) {
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
  Grid mygrad1(gd, n*VectorXd::Ones(gd.NxNyNz));
  Grid mygrad2(gd, n*VectorXd::Ones(gd.NxNyNz));

  printf("Working on f1 of double...\n");
  const double Edouble1 = f1(n);
  printf("Working on f1 of grid...\n");
  const double Egrid1 = f1.integral(nr)/gd.Lat.volume();
  printf("Working on f2 of double...\n");
  const double Edouble2 = f2(n);
  printf("Working on f2 of grid...\n");
  const double Egrid2 = f2.integral(nr)/gd.Lat.volume();

  printf("Edouble1 = %g\n", Edouble1);
  printf("Egrid1   = %g\n", Egrid1);
  printf("Edouble2 = %g\n", Edouble2);
  printf("Egrid2   = %g\n", Egrid2);
  printf("fractional error = %g\n", (Edouble1 - Egrid1)/fabs(Edouble1));
  int retval = 0;
  if (fabs((Edouble1 - Egrid1)/Edouble1) > fraccuracy) {
    printf("FAIL: Error in the energy is too big!\n");
    retval++;
  }
  printf("fractional error3 = %g\n", (Edouble2 - Edouble1)/fabs(Edouble1));
  if (fabs((Edouble2 - Edouble1)/Edouble1) > fraccuracy) {
    printf("FAIL: Error in the energy is too big!\n");
    retval++;
  }
  printf("fractional error2 = %g\n", (Edouble2 - Egrid2)/fabs(Edouble2));
  if (fabs((Edouble2 - Egrid2)/Edouble2) > fraccuracy) {
    printf("FAIL: Error in the energy is too big!\n");
    retval++;
  }

  retval += f1.run_finite_difference_test("f1", nr);
  retval += f2.run_finite_difference_test("f2", nr);
  mygrad1.setZero();
  mygrad2.setZero();
  f1.integralgrad(gd, nr, &mygrad1);
  f2.integralgrad(gd, nr, &mygrad2);
  printf("fractional error in grad = %g\n", (mygrad1[0] - mygrad2[0])/fabs(mygrad2[0]));
  if (fabs((mygrad1[0] - mygrad2[0])/mygrad2[0]) > fraccuracy) {
    printf("FAIL: Error in the grad is too big!\n");
    retval++;
  }
  return retval;
}

int main(int, char **argv) {
  int retval = 0;
  const double kT = 1e-3;
  const FieldFunctional n = EffectivePotentialToDensity(kT);

  FieldFunctional x = Gaussian(1);
  retval += test_functionals("Square vs mul", Pow(2)(x), x*x, 0.1, 1e-13);
  retval += test_functionals("Cube vs mul", Pow(3)(x), x*x*x, 0.1, 1e-12);
  retval += test_functionals("Cube vs other mul", Pow(3)(x), x*(x*x), 0.1, 1e-12);

  // The following tests the chain rule...
  retval += test_functionals("Pow(2)(Pow(2))", Pow(2)(Pow(2)(x)), x*x*x*x, 0.1, 1e-12);

  // The following tests scalar multiply rule...
  retval += test_functionals("Fourth power", 10*Pow(2)(Pow(2)(x)), 10*x*x*x*x, 0.1, 1e-12);

  // The following tests scalar subtraction rule...
  retval += test_functionals("Subtraction",
                             0 - Pow(2)(Pow(2)(x)), -1*Pow(2)(Pow(2)(x)), 0.1, 1e-12);
  retval += test_functionals("Simple subtraction",
                             0 - x, -1*x, 0.1, 1e-12);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
