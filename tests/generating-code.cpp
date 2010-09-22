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
  Functional kT(1e-2), x(Identity()), One(1);
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

  Functional four_pi_r2 = 4*M_PI*sqr(R);
  Functional n2 = ShellConvolve(Rval);
  Functional n3 = StepConvolve(Rval);
  Functional one_minus_n3 = 1 - n3;
  Functional phi1 = (Functional(-1)/four_pi_r2)*n2*log(one_minus_n3);
  phi1.set_name("phi1").create_header("tests/generated/phi1.h", "Phi1", "kT", "R");

  const Functional four_pi_r = 4*M_PI*R;
  Functional n2x = xShellConvolve(Rval);
  Functional n2y = yShellConvolve(Rval);
  Functional n2z = zShellConvolve(Rval);
  Functional phi2 = (sqr(n2) - sqr(n2x) - sqr(n2y) - sqr(n2z))/(four_pi_r*one_minus_n3);
  phi2.set_name("phi2").create_header("tests/generated/phi2.h", "Phi2", "kT", "R");

  Functional phi3rf = n2*(sqr(n2) - 3*(sqr(n2x) + sqr(n2y) + sqr(n2z)))/(24*M_PI*sqr(one_minus_n3));
  phi3rf.set_name("phi3rf").create_header("tests/generated/phi3rf.h", "Phi3rf", "kT", "R");

  //(phi1+phi2+phi3rf).create_header("tests/generated/almostrf.h", "AlmostRF", "kT", "R");

  ((kT*phi1).set_name("phi1")+(kT*phi2).set_name("phi2")+(kT*phi3rf).set_name("phi3")).create_header("tests/generated/almostrf.h", "AlmostRF", "kT", "R");
}
