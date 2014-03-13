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
#include "new/QuadraticFast.h"
#include "new/QuadraticN0Fast.h"
#include "handymath.h"

int check_functional_value(const char *name,
                            const NewFunctional &f,
                            double energy,
                            double fraccuracy = 1e-15) {
  int errors = 0;
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  printf("findiing energy\n");
  double fv = f.energy();
  printf("found energy\n");
  print_double("Energy of Vector:  ", fv);
  printf("\n");
  f.printme("  ");

  if (!(fabs(fv/energy - 1) < fraccuracy)) {
    printf("fv = %.15g\n", fv);
    printf("expected = %.15g\n", energy);
    printf("FAIL: Error in f(n) is %g\n", fv/energy - 1);
    errors++;
  }
  return errors;
}

int main(int, char **argv) {
  int retval = 0;

  int Nx = 100;
  double a = 2;

  Quadratic q(Nx, Nx, Nx);
  q.a1() = a;
  q.a2() = a;
  q.a3() = a;
  q.konst() = 1;
  q.meanval() = 1;
  q.x() = 1.5;
  double energy = 0.5*q.konst()*a*a*a*sqr(q.x()[0]-q.meanval());
  retval += check_functional_value("Quadratic", q, energy, 2e-11);
  retval += q.run_finite_difference_test("Quadratic");

  QuadraticN0 q0(Nx, Nx, Nx);
  q0.a1() = a;
  q0.R() = 1;
  q0.a2() = a;
  q0.a3() = a;
  q0.konst() = 1;
  q0.meanval() = 1;
  q0.x() = 1.5;
  energy = 0.5*q0.konst()*a*a*a*sqr(q0.x()[0]-q0.meanval());
  retval += check_functional_value("Quadratic n0", q0, energy, 2e-11);
  retval += q0.run_finite_difference_test("Quadratic n0");

  for (int i=0; i<Nx*Nx*Nx/2; i++) q0.x()[i] *= 0;
  retval += q0.run_finite_difference_test("Quadratic n0 inhomogeneous");

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
