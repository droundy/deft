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

#include "OptimizedFunctionals.h"
#include "ContactDensity.h"
#include "equation-of-state.h"
#include <time.h>

void took(const char *action) {
  static clock_t start = 0;
  clock_t end = clock();
  printf("    %s took %g seconds.\n", action, (end - double(start))/CLOCKS_PER_SEC);
  start = end;
}

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
  took("Evaluating f1 of n");
  double f2n = f2.integral(water_prop.kT, n);
  took("Evaluating f2 of n");
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
  took("Evaluating grad f1");
  f2.integralgrad(water_prop.kT, n, &gr2);
  took("Evaluating grad f2");
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
  if (fabs(1 - fabs(f1x/f2x)) > fraccuracy) {
    printf("FAIL: Error in double %s is %g as a fraction of %g\n", f1.get_name().c_str(),
           1 - fabs(f1x/f2x), f2x);
    printf("Ratio f1/f2 is %g\n", f1x/f2x);
    errors++;
  } else {
    printf("Functions of double agree!\n");
  }
  
  double f1p = f1.derive(1.2*water_prop.kT, x);
  double f2p = f2.derive(1.2*water_prop.kT, x);
  if (fabs(1 - fabs(f1p/f2p)) > fraccuracy) {
    printf("FAIL: Error in derive double %s is %g as a fraction of %g\n", f1.get_name().c_str(),
           1 - fabs(f1p/f2p), f2p);
    printf("Ratio f1/f2 is %g\n", f1p/f2p);
    errors++;
  } else {
    printf("Derivatives of double agree!\n");
  }

  {
    // The following is redundant, and overkill, but was added when I
    // was desperate to find a bug...

    // Now, let's test that the Functional agrees between its 3D
    // implementation and the homogeneous limit.
    Lattice lat(Cartesian(4.2,0,-0.1), Cartesian(0.1,5.2,0), Cartesian(0,0,5));
    double resolution = 0.2;
    GridDescription gd(lat, resolution);
    const double n = x;
    Grid nr(gd, n*VectorXd::Ones(gd.NxNyNz));

    const double kT = water_prop.kT;
    const double Edouble = f1(kT, n);
    const double Edouble2 = f2(kT, n);
    took("Evaluating functional of a double");
    printf("Edouble = %g\n", Edouble);
    const double Egrid = f1.integral(kT, nr)/gd.Lat.volume();
    const double Egrid2 = f2.integral(kT, nr)/gd.Lat.volume();
    took("Evaluating functional of a 20x20x20 grid");
    
    printf("Edouble %15s = %.10g\n", f1.get_name().c_str(), Edouble);
    printf("Egrid1  %15s = %.10g\n", f1.get_name().c_str(), Egrid);
    printf("Edouble %15s = %.10g\n", f2.get_name().c_str(), Edouble2);
    printf("Egrid   %15s = %.10g\n", f2.get_name().c_str(), Egrid2);
    
    const double deriv_double = f1.derive(kT, n);
    took("Evaluating derivative of a double");
    Grid grad(nr);
    grad.setZero();
    f1.integralgrad(kT, gd, nr, &grad);
    took("Evaluating gradient of a 20x20x20 grid");
    const double deriv_grid = grad[0]/gd.dvolume;
    printf("deriv double = %g\n", deriv_double);
    printf("deriv grid   = %g\n", deriv_grid);
    if (Edouble == 0) {
      // If the true answer is zero, the Edouble should say so!
      printf("absolute error = %g\n", Edouble - Egrid);
      if (fabs(Edouble - Egrid) > fraccuracy) {
        printf("FAIL: Error in the energy is too big!\n");
        errors++;
      }
    } else {
      printf("fractional error in energy = %g\n", (Edouble - Egrid)/fabs(Edouble));
      if (!(fabs((Edouble - Egrid)/Edouble) < fraccuracy)) {
        printf("FAIL: Error in the energy is too big!\n");
        errors++;
      }
    }
    // Check that the derivative agrees...
    if (deriv_double == 0) {
      printf("absolute error = %g\n", deriv_double - deriv_grid);
      if (fabs(deriv_double - deriv_grid) > fraccuracy) {
        printf("FAIL: Error in the gradient is too big!\n");
        errors++;
      }
    } else {
      printf("fractional error in gradient = %g\n", (deriv_double - deriv_grid)/fabs(deriv_double));
      if (!(fabs((deriv_double - deriv_grid)/deriv_double) < fraccuracy)) {
        printf("FAIL: Error in the gradient is just too darn big!\n");
        errors++;
      }
    }
  }
}

int main(int, char **argv) {
  Functional x(Identity());

  Grid n(gd);
  //n = 0.00001*VectorXd::Ones(gd.NxNyNz);
  n = 0.00001*VectorXd::Ones(gd.NxNyNz) + 0.1*(-10*r2(gd)).cwise().exp();

  compare_functionals(HardSpheresNoTensor2(R), HardSpheresNoTensor(R), n, 4e-13);

  compare_functionals(HardSpheresFast(R), HardSpheres(R), n, 1e-13);

  compare_functionals(HardSpheresRFFast(R), HardSpheresRF(R), n, 1e-13);

  compare_functionals(HardSpheresTarazonaFast(R), HardSpheresTarazona(R), n, 1e-13);

  compare_functionals(HardSpheresWBnotensor(R), HardSpheresNoTensor(R), n, 2e-13);

  //compare_functionals(ContactAtSphere(R), ContactDensitySphere(R), n, 3e-14);

  compare_functionals(YuWuCorrelationFast(R), YuWuCorrelation(R), n, 3e-13);

  const double mu = 1e-5;
  Functional f = HardSpheresRFFast(R) + IdealGas() + ChemicalPotential(mu);
  compare_functionals(HardSphereGasRF(R, mu), f, Grid(gd, -water_prop.kT*n.cwise().log()), 1e-12);
 
  f = HardSpheresFast(R) + IdealGas() + ChemicalPotential(mu);
  compare_functionals(HardSphereGas(R, mu), f, Grid(gd, -water_prop.kT*n.cwise().log()), 1e-12);

  double eps = water_prop.epsilonAB;
  double kappa = water_prop.kappaAB;
  double epsdis = 1e-5;
  double lambda = 1.8;
  double lscale = 0.7;
  compare_functionals(SaftFluid(R, eps, kappa, epsdis, lambda, lscale, mu),
                      SaftFluidSlow(R, eps, kappa, epsdis, lambda, lscale, mu),
                      n, 3e-12);

  compare_functionals(SaftFluid2(R, eps, kappa, epsdis, lambda, lscale, mu),
                      SaftFluid(R, eps, kappa, epsdis, lambda, lscale, mu),
                      n, 3e-12);

  //compare_functionals(Association2(R, eps, kappa, epsdis, lambda, lscale),
  //                    Association(R, eps, kappa, epsdis, lambda, lscale),
  //                    n, 3e-12);

  //compare_functionals(Dispersion2(R, epsdis, lambda, lscale),
  //                    Dispersion(R, epsdis, lambda, lscale), n, 3e-12);

  if (errors == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], errors);
  return errors;
}
