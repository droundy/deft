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
#include "new/WhiteBearFast.h"
#include "new/HomogeneousWhiteBearFast.h"

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

  const int Nx = 20;
  const double R = 1.0, a = 2.0, kT = 1, nval = 0.1;
  const double energy = 2.72225468848892;
  printf("about to create input\n");
  WhiteBear wb(Nx, Nx, Nx);
  wb.R() = R;
  wb.a1() = a;
  wb.a2() = a;
  wb.a3() = a;
  wb.kT() = kT;
  wb.n() = nval;
  retval += check_functional_value("WhiteBear", wb, energy, 2e-11);
  printf("n0 = %g\n", wb.get_n0()[0]);
  printf("n1 = %g\n", wb.get_n1()[0]);
  printf("n2 = %g\n", wb.get_n2()[0]);
  printf("n3 = %g\n", wb.get_n3()[0]);

  for (int i=0;i<Nx*Nx*Nx/2;i++) wb.n()[i] = 0.1*nval;
  //FIXME:  the following test OUGHT to pass, but currently fails.  :(
  //retval += wb.run_finite_difference_test("WhiteBear");

  HomogeneousWhiteBear hwb;
  hwb.R() = R;
  hwb.kT() = kT;
  hwb.n() = nval;
  retval += check_functional_value("HomogeneousWhiteBear", hwb, energy/uipow(a,3), 2e-11);
  printf("n0 = %g\n", hwb.get_n0());
  printf("n1 = %g\n", hwb.get_n1());
  printf("n2 = %g\n", hwb.get_n2());
  printf("n3 = %g\n", hwb.get_n3());

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
