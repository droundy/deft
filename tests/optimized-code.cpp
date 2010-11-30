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

int errors = 0;

const double kT = water_prop.kT; // room temperature in Hartree
const double R = 2.7;

double a = 5;
Lattice lat(Cartesian(0,a,a), Cartesian(a,0,a), Cartesian(a,a,0));
//Lattice lat(Cartesian(1.4*rmax,0,0), Cartesian(0,1.4*rmax,0), Cartesian(0,0,1.4*rmax));
GridDescription gd(lat, 0.2);

void compare_functionals(const Functional &f1, const Functional &f2, const Grid &n, double fraccuracy = 1e-15) {
  printf("\n************");
  for (unsigned i=0;i<strlen(f1.get_name());i++) printf("*");
  printf("\n* Testing %s *\n", f1.get_name());
  for (unsigned i=0;i<strlen(f1.get_name());i++) printf("*");
  printf("************\n\n");

  double f1n = f1.integral(n);
  double f2n = f2.integral(n);
  if (fabs(f1n/f2n - 1) > fraccuracy) {
    printf("E1 = %g\n", f1n);
    printf("E2 = %g\n", f2n);
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

  Grid n(gd);
  n = 0.001*VectorXd::Ones(gd.NxNyNz) + 0.001*(-10*r2(gd)).cwise().exp();

  compare_functionals(HardSpheresFast(R, kT), HardSpheres(R, kT), n, 1e-14);

  compare_functionals(HardSpheresRFFast(R, kT), HardSpheresRF(R, kT), n, 1e-14);

  compare_functionals(HardSpheresTarazonaFast(R, kT), HardSpheresTarazona(R, kT), n, 1e-14);

  compare_functionals(HardSpheresWBnotensor(R, kT), HardSpheresNoTensor(R, kT), n, 1e-14);

  
  Functional nn = EffectivePotentialToDensity(kT);
  Functional f = HardSpheres(R, kT) + IdealGas(kT);
  double mu = -f.derive(water_prop.liquid_density);
  f = HardSpheresRFFast(R, kT)(nn) + IdealGasOfVeff(kT) + ChemicalPotential(mu)(nn);
  compare_functionals(HardSphereGasRF(R, kT, mu), f, Grid(gd, -kT*n.cwise().log()), 4e-13);
 
  mu = -(HardSpheres(R, kT)(nn) + IdealGasOfVeff(kT)).derive(water_prop.liquid_density);
  f = HardSpheresFast(R, kT)(nn) + IdealGasOfVeff(kT) + ChemicalPotential(mu)(nn);
  compare_functionals(HardSphereGas(R, kT, mu), f, Grid(gd, -kT*n.cwise().log()), 1e-12);

  double eps = water_prop.epsilonAB;
  double kappa = water_prop.kappaAB;
  double epsdis = 1e-5;
  double lambda = 1.8;
  compare_functionals(SaftFluid(R, kT, eps, kappa, epsdis, lambda, mu),
                      SaftFluidSlow(R, kT, eps, kappa, epsdis, lambda, mu),
                      Grid(gd, -kT*n.cwise().log()), 1e-12);

  if (errors == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], errors);
  return errors;
}
