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
#include "generated/sum.h"
#include "generated/log.h"
#include "generated/log-and-sqr.h"
#include "generated/log-and-sqr-and-inverse.h"
#include "generated/log-one-minus-x.h"
#include "generated/log-one-minus-nbar.h"
#include "generated/sqr-xshell.h"
#include "generated/phi1.h"
#include "generated/phi2.h"
#include "generated/phi3rf.h"
#include "generated/almostrf.h"

int errors = 0;

const double kT = water_prop.kT; // room temperature in Hartree
const double R = 2.7;

double a = 5;
Lattice lat(Cartesian(0,a,a), Cartesian(a,0,a), Cartesian(a,a,0));
//Lattice lat(Cartesian(1.4*rmax,0,0), Cartesian(0,1.4*rmax,0), Cartesian(0,0,1.4*rmax));
GridDescription gd(lat, 0.2);

void compare_functionals(const Functional &f1, const Functional &f2, double fraccuracy = 1e-15) {
  printf("\n************");
  for (unsigned i=0;i<strlen(f1.get_name());i++) printf("*");
  printf("\n* Testing %s *\n", f1.get_name());
  for (unsigned i=0;i<strlen(f1.get_name());i++) printf("*");
  printf("************\n\n");
  
  Grid n(gd);
  n = 0.001*VectorXd::Ones(gd.NxNyNz) + 0.001*(-10*r2(gd)).cwise().exp();

  double f1n = f1.integral(n);
  double f2n = f2.integral(n);
  if (fabs(f1n/f2n - 1) > fraccuracy) {
    printf("FAIL: Error in f(n) is %g\n", f1n/f2n - 1);
    errors++;
  }
  Grid gr1(gd), gr2(gd);
  gr1.setZero();
  gr2.setZero();
  f1.integralgrad(n, &gr1);
  f2.integralgrad(n, &gr2);
  double err = (gr1-gr2).cwise().abs().maxCoeff();
  double mag = gr1.cwise().abs().maxCoeff();
  if (err/mag > fraccuracy) {
    printf("FAIL: Error in grad %s is %g as a fraction of %g\n", f1.get_name(), err/mag, mag);
    errors++;
  }
  errors += f1.run_finite_difference_test(f1.get_name(), n);
  //errors += f2.run_finite_difference_test("other version", n);
}

int main(int, char **argv) {
  Functional x(Identity());
  compare_functionals(Sum(kT), x + kT);

  compare_functionals(Log(), log(x));

  compare_functionals(LogAndSqr(), log(x) + sqr(x));

  compare_functionals(LogAndSqrAndInverse(), log(x) + (sqr(x)-Pow(3)) + Functional(1)/x);

  compare_functionals(LogOneMinusX(), log(1-x));

  compare_functionals(LogOneMinusNbar(R), log(1-StepConvolve(R)));

  compare_functionals(SquareXshell(R), sqr(xShellConvolve(R)));

  const double four_pi_r2 = 4*M_PI*R*R;
  Functional n2 = ShellConvolve(R);
  Functional n3 = StepConvolve(R);
  Functional one_minus_n3 = 1 - n3;
  Functional phi1 = (-1/four_pi_r2)*n2*log(one_minus_n3);
  compare_functionals(Phi1(kT,R), phi1);

  const double four_pi_r = 4*M_PI*R;
  Functional n2x = xShellConvolve(R);
  Functional n2y = yShellConvolve(R);
  Functional n2z = zShellConvolve(R);
  Functional phi2 = (sqr(n2) - sqr(n2x) - sqr(n2y) - sqr(n2z))/(four_pi_r*one_minus_n3);
  compare_functionals(Phi2(kT,R), phi2);

  Functional phi3rf = n2*(sqr(n2) - 3*(sqr(n2x) + sqr(n2y) + sqr(n2z)))/(24*M_PI*sqr(one_minus_n3));
  compare_functionals(Phi3rf(kT,R), phi3rf, 2e-15);

  compare_functionals(AlmostRF(kT,R), kT*(phi1 + phi2 + phi3rf), 2e-15);

  compare_functionals(HardSpheresFast(R, kT), HardSpheres(R, kT), 3e-15);

  compare_functionals(HardSpheresRFFast(R, kT), HardSpheresRF(R,kT), 2e-15);

  compare_functionals(HardSpheresTarazonaFast(R, kT), HardSpheresTarazona(R,kT), 3e-15);

  compare_functionals(HardSpheresWBnotensor(R, kT), HardSpheresNoTensor(R,kT), 2e-15);

  if (errors == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], errors);
  return errors;
}
