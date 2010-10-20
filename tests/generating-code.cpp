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

int main(int, char **) {
  double kTval = 1e-2;
  Functional kT(kTval), x(Identity()), One(1);
  kT.set_name("kT");
  const double Rval = 4;
  Functional R(Rval);
  R.set_name("R");
  (kT + x).create_header("tests/generated/sum.h", "Sum", "kT");
  log(x).create_header("tests/generated/log.h", "Log");
  (log(x)+sqr(x)).create_header("tests/generated/log-and-sqr.h", "LogAndSqr");
  (One + sqr(x) - Pow(3) + One/x).create_header("tests/generated/one-over-x.h", "OneOverX");
  (log(x)+(sqr(x)-Pow(3))+One/x).create_header("tests/generated/log-and-sqr-and-inverse.h", "LogAndSqrAndInverse");
  log(1-x).create_header("tests/generated/log-one-minus-x.h", "LogOneMinusX");
  log(1-StepConvolve(Rval)).create_header("tests/generated/log-one-minus-nbar.h", "LogOneMinusNbar", "R");
  sqr(xShellConvolve(Rval)).create_header("tests/generated/sqr-xshell.h", "SquareXshell", "R");

  Functional veff = EffectivePotentialToDensity(kTval);

  sqr(veff).create_header("tests/generated/sqr-Veff.h", "SquareVeff", "kT", "R");

  IdealGasOfVeff(kTval).create_header("tests/generated/ideal-gas.h", "IdealGasFast", "kT");

  Functional n2 = ShellConvolve(Rval);
  Functional n3 = StepConvolve(Rval);
  (sqr(n2)+sqr(n3)).set_name("n2_and_n3").create_header("tests/generated/n2_and_n3.h", "n2_and_n3", "R");

  Functional four_pi_r2 = 4*M_PI*sqr(R);
  Functional one_minus_n3 = 1 - n3;
  Functional phi1 = (Functional(-1)/four_pi_r2)*n2*log(one_minus_n3);
  phi1.set_name("phi1").create_header("tests/generated/phi1.h", "Phi1", "kT", "R");

  phi1(veff).create_header("tests/generated/phi1-Veff.h", "Phi1Veff", "kT", "R");

  const Functional four_pi_r = 4*M_PI*R;
  Functional n2x = xShellConvolve(Rval);
  Functional n2y = yShellConvolve(Rval);
  Functional n2z = zShellConvolve(Rval);
  Functional phi2 = (sqr(n2) - sqr(n2x) - sqr(n2y) - sqr(n2z))/(four_pi_r*one_minus_n3);
  phi2.set_name("phi2").create_header("tests/generated/phi2.h", "Phi2", "kT", "R");

  phi2(veff).create_header("tests/generated/phi2-Veff.h", "Phi2Veff", "kT", "R");

  Functional phi3rf = n2*(sqr(n2) - 3*(sqr(n2x) + sqr(n2y) + sqr(n2z)))/(24*M_PI*sqr(one_minus_n3));
  phi3rf.set_name("phi3rf").create_header("tests/generated/phi3rf.h", "Phi3rf", "kT", "R");

  phi3rf(veff).create_header("tests/generated/phi3rf-Veff.h", "Phi3rfVeff", "kT", "R");

  (phi1(veff) + IdealGasOfVeff(kTval) + ChemicalPotential(0)(veff)).create_header("tests/generated/phi1-plus.h",
                                                                                  "Phi1plus", "R", "kT", "mu");

  //(phi1+phi2+phi3rf).create_header("tests/generated/almostrf.h", "AlmostRF", "kT", "R");

  ((kT*phi1).set_name("phi1")+(kT*phi2).set_name("phi2")+(kT*phi3rf).set_name("phi3")).create_header("tests/generated/almostrf.h", "AlmostRF", "kT", "R");

  (phi1+phi2+phi3rf).create_header("tests/generated/almostrfnokt.h", "AlmostRFnokT", "kT", "R");
}
