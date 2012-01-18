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
#include "equation-of-state.h"

#include "generated-haskell/nice-sum.h"
#include "generated-haskell/nice-quad.h"
#include "generated-haskell/nice-sqrt.h"
#include "generated-haskell/nice-sqrtandmore.h"
#include "generated-haskell/nice-log.h"
#include "generated-haskell/nice-logandsqr.h"
#include "generated-haskell/nice-logandsqrandinverse.h"
#include "generated-haskell/nice-logoneminusx.h"
#include "generated-haskell/nice-nbar.h"

#include "generated/sum.h"
#include "generated/quadratic.h"
#include "generated/sqrt.h"
#include "generated/sqrt-and-more.h"
#include "generated/log.h"
#include "generated/log-and-sqr.h"
#include "generated/log-and-sqr-and-inverse.h"
#include "generated/log-one-minus-x.h"
#include "generated/nbar.h"
#include "generated/log-one-minus-nbar.h"
#include "generated/sqr-xshell.h"
#include "generated/n2_and_n3.h"
#include "generated/sqr-Veff.h"
#include "generated/ideal-gas.h"

#include "generated/phi1.h"
#include "generated/phi1-Veff.h"
#include "generated/phi1-plus.h"
#include "generated/phi2.h"
#include "generated/phi2-Veff.h"
#include "generated/phi3rf.h"
#include "generated/phi3rf-Veff.h"
#include "generated/almostrf.h"
#include "generated/almostrfnokt.h"

#include "handymath.h"

int errors = 0;

double a = 5;
Lattice lat(Cartesian(0,a,a), Cartesian(a,0,a), Cartesian(a,a,0));
//Lattice lat(Cartesian(1.4*rmax,0,0), Cartesian(0,1.4*rmax,0), Cartesian(0,0,1.4*rmax));
GridDescription gd(lat, 0.2);

void compare_functionals(const Functional &f1, const Functional &f2,
                         double kT, const Grid &n, double fraccuracy = 1e-15,
                         double x = 0.001, double fraccuracydoub = 1e-15) {
  printf("\n************");
  for (unsigned i=0;i<f1.get_name().size();i++) printf("*");
  printf("\n* Testing %s *\n", f1.get_name().c_str());
  for (unsigned i=0;i<f1.get_name().size();i++) printf("*");
  printf("************\n\n");

  printf("First energy:\n");
  double f1n = f1.integral(kT, n);
  print_double("first energy is:               ", f1n);
  printf("\n");
  f1.print_summary("", f1n, f1.get_name().c_str());
  printf("Second energy:\n");
  double f2n = f2.integral(kT, n);
  print_double("second energy is:              ", f2n);
  printf("\n");
  f2.print_summary("", f2n, f2.get_name().c_str());
  if (fabs(f1n/f2n - 1) > fraccuracy) {
    printf("E1 = %g\n", f1n);
    printf("E2 = %g\n", f2n);
    printf("FAIL: Error in f(n) is %g\n", f1n/f2n - 1);
    errors++;
  }
  Grid gr1(gd), gr2(gd);
  gr1.setZero();
  gr2.setZero();
  f1.integralgrad(kT, n, &gr1);
  f2.integralgrad(kT, n, &gr2);
  double err = (gr1-gr2).cwise().abs().maxCoeff();
  double mag = gr1.cwise().abs().maxCoeff();
  if (err/mag > fraccuracy) {
    printf("FAIL: Error in grad %s is %g as a fraction of %g\n", f1.get_name().c_str(), err/mag, mag);
    errors++;
  }
  errors += f1.run_finite_difference_test(f1.get_name().c_str(), kT, n);

  double f1x = f1(kT, x);
  double f2x = f2(kT, x);
  if (1 - fabs(f1x/f2x) > fraccuracydoub) {
    printf("FAIL: Error in double %s is %g as a fraction of %g it is %g\n", f1.get_name().c_str(),
           1 - fabs(f1x/f2x), f2x, f1x);
    errors++;
  }
  
  double f1p = f1.derive(kT, x);
  double f2p = f2.derive(kT, x);
  if (1 - fabs(f1p/f2p) > fraccuracydoub) {
    printf("FAIL: Error in derive double %s is %g as a fraction of %g\n", f1.get_name().c_str(),
           1 - fabs(f1p/f2p), f2p);
    errors++;
  }
  
  //errors += f2.run_finite_difference_test("other version", n);
}

int main(int, char **argv) {
  const double myT = water_prop.kT; // room temperature in Hartree
  const double R = 2.7;

  Functional x(Identity());
  Grid n(gd);
  n = 0.001*VectorXd::Ones(gd.NxNyNz) + 0.001*(-10*r2(gd)).cwise().exp();

  compare_functionals(NiceSum(), x + kT, myT, n, 4e-13);

  compare_functionals(NiceQuad(), sqr(x + kT) - x + 2*kT, myT, n, 2e-12);

  compare_functionals(NiceSqrt(), sqrt(x) , myT, n, 1e-12);

  compare_functionals(NiceSqrtandMore(), sqrt(x) - x + 2*kT, myT, n, 1e-12);

  compare_functionals(NiceLog(), log(x), myT, n, 7e-14);

  compare_functionals(NiceLogandSqr(), log(x) + sqr(x), myT, n, 6e-13);

  compare_functionals(NiceLogandSqrandInverse(), log(x) + (sqr(x)-Pow(3)) + Functional(1)/x, myT, n, 3e-10);

  compare_functionals(NiceLogOneMinusX(), log(1-x), myT, n, 1e-12);

  compare_functionals(NiceNbar(R), StepConvolve(R), myT, n, 3e-13);

  compare_functionals(Sum(), x + kT, myT, n, 2e-13);

  compare_functionals(Quadratic(), sqr(x + kT) - x + 2*kT, myT, n, 2e-12);

  compare_functionals(Sqrt(), sqrt(x), myT, n, 1e-12);

  compare_functionals(SqrtAndMore(), sqrt(x) - x + 2*kT, myT, n, 1e-12);

  compare_functionals(LogAndSqr(), log(x) + sqr(x), myT, n, 3e-14);

  compare_functionals(LogAndSqrAndInverse(), log(x) + (sqr(x)-Pow(3)) + Functional(1)/x, myT, n, 3e-10);

  compare_functionals(LogOneMinusX(), log(1-x), myT, n, 1e-12);

  compare_functionals(LogOneMinusNbar(R), log(1-StepConvolve(R)), myT, n, 1e-13);

  compare_functionals(Nbar(R), StepConvolve(R), myT, n, 3e-13);

  compare_functionals(SquareXshell(R), sqr(xShellConvolve(R)), myT, n);

  Functional n2 = ShellConvolve(R);
  Functional n3 = StepConvolve(R);
  compare_functionals(n2_and_n3(R), sqr(n2) + sqr(n3), myT, n, 1e-14);

  const double four_pi_r2 = 4*M_PI*R*R;
  Functional one_minus_n3 = 1 - n3;
  Functional phi1 = (-1/four_pi_r2)*n2*log(one_minus_n3);
  compare_functionals(Phi1(R), phi1, myT, n, 1e-13);

  const double four_pi_r = 4*M_PI*R;
  Functional n2x = xShellConvolve(R);
  Functional n2y = yShellConvolve(R);
  Functional n2z = zShellConvolve(R);
  Functional phi2 = (sqr(n2) - sqr(n2x) - sqr(n2y) - sqr(n2z))/(four_pi_r*one_minus_n3);
  compare_functionals(Phi2(R), phi2, myT, n, 1e-14);

  Functional phi3rf = n2*(sqr(n2) - 3*(sqr(n2x) + sqr(n2y) + sqr(n2z)))/(24*M_PI*sqr(one_minus_n3));
  compare_functionals(Phi3rf(R), phi3rf, myT, n, 1e-13);

  compare_functionals(AlmostRF(R), myT*(phi1 + phi2 + phi3rf), myT, n, 2e-14);

  Functional veff = EffectivePotentialToDensity();
  compare_functionals(SquareVeff(R), sqr(veff), myT, Grid(gd, -myT*n.cwise().log()), 1e-12);

  compare_functionals(AlmostRFnokT(R), phi1 + phi2 + phi3rf, myT, n, 3e-14);

  compare_functionals(AlmostRF(R),
                      (myT*phi1).set_name("phi1") + (myT*phi2).set_name("phi2") + (myT*phi3rf).set_name("phi3"),
                      myT, n, 4e-14);

  compare_functionals(Phi1Veff(R), phi1(veff), myT, Grid(gd, -myT*n.cwise().log()), 1e-13);

  compare_functionals(Phi2Veff(R), phi2(veff), myT, Grid(gd, -myT*n.cwise().log()), 1e-14);

  compare_functionals(Phi3rfVeff(R), phi3rf(veff), myT, Grid(gd, -myT*n.cwise().log()), 1e-13, 0.001, 2e-15);

  compare_functionals(IdealGasFast(), IdealGasOfVeff, myT, Grid(gd, -myT*n.cwise().log()),
                      1e-12);

  double mu = -1;
  compare_functionals(Phi1plus(R, mu),
                      phi1(veff) + IdealGasOfVeff + ChemicalPotential(mu)(veff),
                      myT, Grid(gd, -myT*n.cwise().log()), 1e-12);

  if (errors == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], errors);
  return errors;
}
