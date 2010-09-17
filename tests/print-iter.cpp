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

Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,1));
GridDescription gd(lat, 2, 2, 2);

int test_print(const char *name, Functional f) {
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  int retval = 0;

  printf("Energy is %g\n", f.integral(gd, VectorXd::Ones(gd.NxNyNz)));
  f.print_iteration("\tPREFIX:", 1);

  return retval;
}

int main(int, char **argv) {
  int retval = 0;

  Functional x(Identity());

  retval += test_print("sqr", sqr(x).set_name("foobar"));

  Functional sqrandlog = sqr(x).set_name("sqr") + log(x).set_name("log");
  Functional sqronly = sqr(x).set_name("sqronly");
  retval += test_print("sqr + log", sqrandlog);
  retval += test_print("(sqr + log)(sqr)", sqrandlog(sqronly));

  retval += test_print("sqr(sqr(x)) + sqr", sqr(sqr(x)).set_name("foobar") + sqr(x).set_name("bazbar"));

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
