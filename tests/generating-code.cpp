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

#include "Functionals.h"
#include "utilities.h"
#include "CallMe.h"

int main(int, char **) {
  Functional x(Identity()), One(1);
  const double Rval = 4;
  Functional R(Rval, "R");
  const char *noargs[] = { 0, 0 };
  (kT.set_name("temp") + x).create_header("tests/generated/sum.h", "Sum", noargs);
  (sqr(kT + x) - x + 2*kT).create_header("tests/generated/quadratic.h", "Quadratic", noargs);
  sqrt(x).create_header("tests/generated/sqrt.h", "Sqrt", noargs);
  (CallMe(sqrt(x), "Sqrt", "()") - x + 2*kT).create_header("tests/generated/sqrt-and-more.h",
                                                           "SqrtAndMore", noargs);
  log(x).create_header("tests/generated/log.h", "Log", noargs);
  (log(x)+sqr(x)).create_header("tests/generated/log-and-sqr.h", "LogAndSqr", noargs);
  (One + sqr(x) - Pow(3) + One/x).create_header("tests/generated/one-over-x.h", "OneOverX", noargs);
  (log(x)+(sqr(x)-Pow(3))+One/x).create_header("tests/generated/log-and-sqr-and-inverse.h",
                                               "LogAndSqrAndInverse", noargs);
  log(1-x).create_header("tests/generated/log-one-minus-x.h", "LogOneMinusX", noargs);
  const char *R_arg[] = { "R", 0 };
  const char *R_mu_arg[] = { "R", "mu", 0 };
  StepConvolve(Rval).create_header("tests/generated/nbar.h", "Nbar", R_arg);
  log(1-StepConvolve(Rval)).create_header("tests/generated/log-one-minus-nbar.h",
                                          "LogOneMinusNbar", R_arg);
  sqr(xShellConvolve(Rval)).create_header("tests/generated/sqr-xshell.h", "SquareXshell", R_arg);

  Functional veff = EffectivePotentialToDensity();

  sqr(veff).create_header("tests/generated/sqr-Veff.h", "SquareVeff", R_arg);

  IdealGasOfVeff().create_header("tests/generated/ideal-gas.h", "IdealGasFast", noargs);

  Functional n2 = ShellConvolve(Rval);
  Functional n3 = StepConvolve(Rval);
  (sqr(n2)+sqr(n3)).set_name("n2_and_n3").create_header("tests/generated/n2_and_n3.h", "n2_and_n3", R_arg);

  Functional four_pi_r2 = 4*M_PI*sqr(R);
  Functional one_minus_n3 = 1 - n3;
  Functional phi1 = (Functional(-1)/four_pi_r2)*n2*log(one_minus_n3);
  //phi1 = n2*log(one_minus_n3);
  //phi1 = (Functional(-1)/four_pi_r2)*n2;
  //phi1 = (Functional(-1)/four_pi_r2);
  phi1.set_name("phi1").create_header("tests/generated/phi1.h", "Phi1", R_arg);

  phi1(veff).create_header("tests/generated/phi1-Veff.h", "Phi1Veff", R_arg);

  const Functional four_pi_r = 4*M_PI*R;
  Functional n2x = xShellConvolve(Rval);
  Functional n2y = yShellConvolve(Rval);
  Functional n2z = zShellConvolve(Rval);
  Functional phi2 = (sqr(n2) - sqr(n2x) - sqr(n2y) - sqr(n2z))/(four_pi_r*one_minus_n3);
  phi2.set_name("phi2").create_header("tests/generated/phi2.h", "Phi2", R_arg);

  phi2(veff).create_header("tests/generated/phi2-Veff.h", "Phi2Veff", R_arg);

  Functional phi3rf = n2*(sqr(n2) - 3*(sqr(n2x) + sqr(n2y) + sqr(n2z)))/(24*M_PI*sqr(one_minus_n3));
  phi3rf.set_name("phi3rf").create_header("tests/generated/phi3rf.h", "Phi3rf", R_arg);

  phi3rf(veff).create_header("tests/generated/phi3rf-Veff.h", "Phi3rfVeff", R_arg);

  (phi1(veff) + IdealGasOfVeff() + ChemicalPotential(0)(veff)).create_header("tests/generated/phi1-plus.h",
                                                                           "Phi1plus", R_mu_arg);

  //(phi1+phi2+phi3rf).create_header("tests/generated/almostrf.h", "AlmostRF", "R");

  ((kT*phi1).set_name("phi1")+(kT*phi2).set_name("phi2")+(kT*phi3rf).set_name("phi3")).create_header("tests/generated/almostrf.h", "AlmostRF", R_arg);

  (CallMe(phi1, "Phi1", "(R)")+phi2+phi3rf).create_header("tests/generated/almostrfnokt.h", "AlmostRFnokT", R_arg);
}
