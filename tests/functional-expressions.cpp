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

  if (e.printme() != expected) {
    printf("FAIL: Got instead: \"%s\"\n", e.printme().c_str());
    retval++;
  }
}


int main(int, char **argv) {
  const double kT = 1e3, R = 2.7;
  const double four_pi_r2 = 4*M_PI*R*R;
  Functional n2 = ShellConvolve(R);
  Functional n3 = StepConvolve(R);
  Functional one_minus_n3 = 1 - n3;
  
  test_expression(IdealGas(kT).printme(Expression("x")),
                  "kT*choose(1e-90, -207.2326583694641*(x - 1e-90*VectorXd::Ones(gd.NxNyNz)) - 2.082326583694641e-88*VectorXd::Ones(gd.NxNyNz), x.cwise()*x.cwise().log() - x)");

  test_expression(sqr(xShellConvolve(R)).printme(Expression("x")),
                  "xShellConvolve(R)(x).cwise().square()");

  test_expression(sqr(xShellConvolve(R)).grad(dV, false).printme(Expression("x")),
                  "2*xShellConvolve(R)(-(xShellConvolve(R)(x)*gd.dvolume))");

  
  test_expression(((-1/four_pi_r2)*n2).printme(Expression("n")),
                  "-0.01091597689244824*ifft(gd, shell(gd, R).cwise()*fft(gd, n))");

  Functional phi1 = (-1/four_pi_r2)*n2*log(one_minus_n3);
  test_expression(phi1.printme(Expression("n")),
                  "(-0.01091597689244824*ifft(gd, shell(gd, R).cwise()*fft(gd, n))).cwise()*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, n))).cwise().log()");

  //  test_expression(HardSpheres(kT, R).printme(Expression("n")),
  //                  "");

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
