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

void test_expression_type(const char *name, const Expression &e, const char *expected) {
  printf("\n**********************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing type of \"%s\" *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("**********************\n");

  if (!e.typeIs(expected)) {
    printf("Expected: \"%s\"\n", expected);
    printf("FAIL: Got instead: \"%s\"\n", e.ctype());
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
  const double R = 2.7;
  Functional RR(R, "R");
  const Functional four_pi_r2 = 4*M_PI*sqr(RR);
  Functional x = Identity();
  Functional n2 = ShellConvolve(R);
  Functional n2x = xShellConvolve(R);
  Functional n3 = StepConvolve(R);
  Functional one_minus_n3 = 1 - n3;
  
  test_expression("sqr(n2)", sqr(n2).printme(Expression("kT"), Expression("x")),
                  "ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise().square()");
  test_expression("Pow(2)(n2)", Pow(2)(n2).printme(Expression("kT"), Expression("x")),
                  "ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise().square()");
  
  test_expression("IdealGasOfVeff",
                  IdealGasOfVeff.printme(Expression("kT"), Expression("x")),
                  "(-(x + kT.cwise()*(((32811.993*kT/6.283185307179586).cwise()*(32811.993*kT/6.283185307179586).cwise().square()).cwise()*(32811.993*kT/6.283185307179586).cwise().sqrt()).cwise().log() + kT)).cwise()*((-x).cwise()/kT).cwise().exp()");

  test_expression("kT*xxx",
                  (kT*sqr(xShellConvolve(R))).printme(Expression("kT").set_type("double"),
                                                      Expression("x")),
                  "kT*ifft(gd, xshell(gd, R).cwise()*fft(gd, x)).cwise().square()");

  test_expression("sqr(n1)",
                  sqr(xShellConvolve(R)).grad(dV, Identity(), false).printme(Expression("kT"), Expression("x")),
                  "-2*gd.dvolume*ifft(gd, xshell(gd, R).cwise()*fft(gd, ifft(gd, xshell(gd, R).cwise()*fft(gd, x))))");


  test_expression("foobar",
                  (four_pi_r2*Identity()).printme(Expression("kT"), Expression("x")), "12.56637061435917*(R*R)*x");
  
  test_expression("foobar",
                  ((Functional(-1)/four_pi_r2)*n2).printme(Expression("kT"), Expression("n")),
                  "-1/(12.56637061435917*(R*R))*ifft(gd, shell(gd, R).cwise()*fft(gd, n))");

  Functional phi1 = (Functional(-1)/four_pi_r2)*n2*log(one_minus_n3);
  test_expression("phi1",
                  phi1.printme(Expression("kT"), Expression("n")),
                  "(-1/(12.56637061435917*(R*R))*ifft(gd, shell(gd, R).cwise()*fft(gd, n))).cwise()*(VectorXd::Ones(gd.NxNyNz) + -ifft(gd, step(gd, R).cwise()*fft(gd, n))).cwise().log()");

  //  test_expression(HardSpheres(kT, R).printme(Expression("kT"), Expression("n")),
  //                  "");

  test_expression_type("phi1", phi1.printme(Expression("kT"), Expression("n")), "Grid");
  test_expression_type("phi1 of double",
                       phi1.printme(Expression("kT").set_type("double"),
                                    Expression("n").set_type("double")), "double");


  Expression n2sqr(sqr(n2).printme(Expression("kT"), Expression("x")));
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

  test_expression("n3(double)", n3.printme(Expression("kT"),
                                           Expression("x").set_type("double")),
                  "4.188790204786391*(R*R*R)*x");
  test_expression("n2(double)", n2.printme(Expression("kT"), Expression("x").set_type("double")),
                  "12.56637061435917*(R*R)*x");
  test_expression("xshell(double)", n2x.printme(Expression("kT"), Expression("x").set_type("double")),
                  "0");

  test_expression_type("exp(x/(-kT))", exp(x/-kT).printme(Expression("kT"), Expression("x")),
                       "Grid");

  test_expression_type("exp(x/(-kT))", exp(x/Functional(-kT)).printme(Expression("kT").set_type("double"),
                                                          Expression("x").set_type("double")),
                       "double");

  test_expression_type("Step(exp(x/(-kT)))",
                       (StepConvolve(R)(exp(x/Functional(-kT)))).printme(Expression("kT").set_type("double"),
                                                             Expression("x").set_type("double")),
                       "double");

  test_expression_type("Step(exp(x/(-kT)))",
                       (StepConvolve(R)(exp(x/Functional(-kT)))).printme(Expression("kT"), Expression("x")),
                       "Grid");

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
