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
#include <string.h>
#include "Functionals.h"

int retval = 0;

void test_expression(const Expression &e, const char *expected) {
  printf("\n**************");
  for (unsigned i=0;i<strlen(expected);i++) printf("*");
  printf("\n* Testing \"%s\" *\n", expected);
  for (unsigned i=0;i<strlen(expected);i++) printf("*");
  printf("**************\n");
  int retval = 0;

  if (e.printme() != expected) {
    printf("FAIL: Got instead: \"%s\"\n", e.printme().c_str());
    retval++;
  }
}


int main(int, char **argv) {

  test_expression(IdealGas(1e-3).printme(Expression("x")),
                  "choose(1e-90, ((x + -1e-90)*-207.233 + -2.08233e-88)*0.001, (x*x.log() - x)*0.001)");

  test_expression(sqr(xShellConvolve(3)).printme(Expression("x")),
                  "VShellConvolve(R)(x).square()");

  test_expression(sqr(xShellConvolve(3)).grad(dV, false).printme(Expression("x")),
                  "VShellConvolve(R)(-VShellConvolve(R)(x)*gd.dvolume)*2");

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
