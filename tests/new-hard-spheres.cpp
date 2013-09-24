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

  if (fabs(fv/energy - 1) > fraccuracy) {
    printf("fv = %g\n", fv);
    printf("expected = %g\n", energy);
    printf("FAIL: Error in f(n) is %g\n", fv/energy - 1);
    errors++;
  }
  return errors;
}

int main(int, char **argv) {
  int retval = 0;

  const int Nx = 100;
  const double R = 1.0, a = 5.0, kT = 1, nval = 0.1;
  const double energy = 42.53522950699669281;
  printf("about to create input\n");
  WhiteBear wb(Nx, Nx, Nx);
  wb.R() = R;
  wb.a1() = a;
  wb.a2() = a;
  wb.a3() = a;
  wb.kT() = kT;
  wb.n() = nval;
  retval += check_functional_value("WhiteBear", wb, energy);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
