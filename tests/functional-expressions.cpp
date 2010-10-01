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

void test_expression(const char *name, const Expression &e, const char *expected) {
  printf("\n**************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing \"%s\" *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("**************\n");

  if (e.printme() != expected) {
    printf("Expected: \"%s\"\n", expected);
    printf("FAIL: Got instead: \"%s\"\n", e.printme().c_str());
    retval++;
  }
}

void test_subexpression(const Expression &e, const char *expected) {
  printf("\n******************");
  for (unsigned i=0;i<strlen(expected);i++) printf("*");
  printf("\n* Testing CSE \"%s\" *\n", expected);
  for (unsigned i=0;i<strlen(expected);i++) printf("*");
  printf("******************\n");

  if (e.FindCommonSubexpression().printme() != expected) {
    printf("FAIL: Got instead: \"%s\"\n", e.FindCommonSubexpression().printme().c_str());
    retval++;
  }
}

int main(int, char **argv) {
  const double kT = 1e-3, R = 2.7;
  Functional RR(R);
  RR.set_name("R");
  const Functional four_pi_r2 = 4*M_PI*sqr(RR);
  Functional x = Identity();
  Functional n2 = ShellConvolve(R);
  Functional n3 = StepConvolve(R);
  Functional one_minus_n3 = 1 - n3;
  
  test_expression("sqr(n2)", sqr(n2).printme(Expression("x")),
                  "ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise().square()");
  test_expression("Pow(2)(n2)", Pow(2)(n2).printme(Expression("x")),
                  "ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise().square()");
  
  test_expression("IdealGas",
                  IdealGas(kT).printme(Expression("x")),
                  "kT*choose(1e-90, -207.2326583694641*(x - 1e-90*VectorXd::Ones(gd.NxNyNz)) - 2.082326583694641e-88*VectorXd::Ones(gd.NxNyNz), x.cwise()*x.cwise().log() - x)");

  test_expression("kT*xxx",
                  (Functional(kT).set_name("kT")*sqr(xShellConvolve(R))).printme(Expression("x")),
                  "kT*ifft(gd, xshell(gd, R).cwise()*fft(gd, x)).cwise().square()");

  test_expression("sqr(n1)",
                  sqr(xShellConvolve(R)).grad(dV, Identity(), false).printme(Expression("x")),
                  "ifft(gd, xshell(gd, R).cwise()*fft(gd, -(2*ifft(gd, xshell(gd, R).cwise()*fft(gd, x))*gd.dvolume)))");

  test_expression("foobar",
                  (four_pi_r2*Identity()).printme(Expression("x")), "12.56637061435917*(R*R)*x");
  
  test_expression("foobar",
                  ((Functional(-1)/four_pi_r2)*n2).printme(Expression("n")),
                  "-1/(12.56637061435917*(R*R))*ifft(gd, shell(gd, R).cwise()*fft(gd, n))");

  Functional phi1 = (Functional(-1)/four_pi_r2)*n2*log(one_minus_n3);
  test_expression("phi1",
                  phi1.printme(Expression("n")),
                  "(-1/(12.56637061435917*(R*R))*ifft(gd, shell(gd, R).cwise()*fft(gd, n))).cwise()*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, n))).cwise().log()");

  //  test_expression(HardSpheres(kT, R).printme(Expression("n")),
  //                  "");

  Expression n2sqr(sqr(n2).printme(Expression("x")));
  test_subexpression(n2sqr, "fft(gd, x)");

  Expression n2sqr_simple(n2sqr);
  n2sqr_simple.EliminateThisSubexpression(n2sqr.FindCommonSubexpression(), "fftx");
  test_expression("n2sqr_simple",
                  n2sqr_simple, "ifft(gd, shell(gd, R).cwise()*fftx).cwise().square()");

  Expression n2sqr_simpler = n2sqr_simple;
  test_subexpression(n2sqr_simpler, "ifft(gd, shell(gd, R).cwise()*fftx)");

  n2sqr_simpler.EliminateThisSubexpression(n2sqr_simpler.FindCommonSubexpression(), "n2");
  test_expression("n2sqr_simpler",
                  n2sqr_simpler, "n2.cwise().square()");
  

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
