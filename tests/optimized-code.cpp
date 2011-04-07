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

int errors = 0;

const double R = 2.7;

double a = 5;
Lattice lat(Cartesian(0,a,a), Cartesian(a,0,a), Cartesian(a,a,0));
//Lattice lat(Cartesian(1.4*rmax,0,0), Cartesian(0,1.4*rmax,0), Cartesian(0,0,1.4*rmax));
GridDescription gd(lat, 0.2);

void compare_functionals(const Functional &f1, const Functional &f2, const Grid &n, double fraccuracy = 1e-15) {
  printf("\n************");
  for (unsigned i=0;i<f1.get_name().size();i++) printf("*");
  printf("\n* Testing %s *\n", f1.get_name().c_str());
  for (unsigned i=0;i<f1.get_name().size();i++) printf("*");
  printf("************\n\n");

  double f1n = f1.integral(water_prop.kT, n);
  double f2n = f2.integral(water_prop.kT, n);
  if (fabs(f1n/f2n - 1) > fraccuracy) {
    printf("E1 = %g\n", f1n);
    printf("E2 = %g\n", f2n);
    printf("FAIL: Error in f(n) is %g\n", f1n/f2n - 1);
    errors++;
  }
  Grid gr1(gd), gr2(gd);
  gr1.setZero();
  gr2.setZero();
  f1.integralgrad(water_prop.kT, n, &gr1);
  f2.integralgrad(water_prop.kT, n, &gr2);
  double err = (gr1-gr2).cwise().abs().maxCoeff();
  double mag = gr1.cwise().abs().maxCoeff();
  if (err/mag > fraccuracy) {
    printf("FAIL: Error in grad %s is %g as a fraction of %g\n", f1.get_name().c_str(), err/mag, mag);
    errors++;
  }
  errors += f1.run_finite_difference_test(f1.get_name().c_str(), water_prop.kT, n);
  //errors += f2.run_finite_difference_test("other version", n);
  
  const double x = 0.001;
  double f1x = f1(1.2*water_prop.kT, x);
  double f2x = f2(1.2*water_prop.kT, x);
  if (1 - fabs(f1x/f2x) > fraccuracy) {
    printf("FAIL: Error in double %s is %g as a fraction of %g\n", f1.get_name().c_str(),
           1 - fabs(f1x/f2x), f2x);
    errors++;
  }
  
  double f1p = f1.derive(1.2*water_prop.kT, x);
  double f2p = f2.derive(1.2*water_prop.kT, x);
  if (1 - fabs(f1p/f2p) > fraccuracy) {
    printf("FAIL: Error in derive double %s is %g as a fraction of %g\n", f1.get_name().c_str(),
           1 - fabs(f1p/f2p), f2p);
    errors++;
  }
}

int main(int, char **argv) {
  Functional x(Identity());

  Grid n(gd);
  n = 0.00001*VectorXd::Ones(gd.NxNyNz) + 0.001*(-10*r2(gd)).cwise().exp();

  compare_functionals(HardSpheresFast(R), HardSpheres(R), n, 2e-14);

  compare_functionals(HardSpheresRFFast(R), HardSpheresRF(R), n, 2e-14);

  compare_functionals(HardSpheresTarazonaFast(R), HardSpheresTarazona(R), n, 2e-14);

  compare_functionals(HardSpheresWBnotensor(R), HardSpheresNoTensor(R), n, 2e-14);

  
  Functional nn = EffectivePotentialToDensity();
  double mu = find_chemical_potential(HardSpheres(R)(nn) + IdealGasOfVeff, water_prop.kT,
                                      water_prop.liquid_density);
  Functional f = HardSpheresRFFast(R)(nn) + IdealGasOfVeff + ChemicalPotential(mu)(nn);
  compare_functionals(HardSphereGasRF(R, mu), f, Grid(gd, -water_prop.kT*n.cwise().log()), 1e-12);
 
  f = HardSpheresFast(R)(nn) + IdealGasOfVeff + ChemicalPotential(mu)(nn);
  compare_functionals(HardSphereGas(R, mu), f, Grid(gd, -water_prop.kT*n.cwise().log()), 1e-12);

  double eps = water_prop.epsilonAB;
  double kappa = water_prop.kappaAB;
  double epsdis = 1e-5;
  double lambda = 1.8;
  double lscale = 0.7;
  compare_functionals(SaftFluid(R, eps, kappa, epsdis, lambda, lscale, mu),
                      SaftFluidSlow(R, eps, kappa, epsdis, lambda, lscale, mu),
                      Grid(gd, -water_prop.kT*n.cwise().log()), 1e-12);

  if (errors == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], errors);
  return errors;
}
