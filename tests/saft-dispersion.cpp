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
#include "Functionals.h"

int retval = 0;

void test_functional(const char *name, Functional f, double n, double etrue,
                     double fraccuracy=1e-14) {
  printf("\n**************************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s of %10g *\n", name, n);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("**************************\n\n");

  // Here we set up the lattice.
  //Lattice lat(Cartesian(4.2,0,-0.1), Cartesian(0.1,5.2,0), Cartesian(0,0,5));
  Lattice lat(Cartesian(1,0,0), Cartesian(0,1,0), Cartesian(0,0,1));
  double resolution = 0.2;
  GridDescription gd(lat, resolution);
  Grid nr(gd, n*VectorXd::Ones(gd.NxNyNz));

  const double Edouble = f(n);
  const double Egrid = f.integral(nr)/gd.Lat.volume();
  f.print_summary("", Egrid);

  printf("Edouble = %g\n", Edouble);
  printf("Egrid   = %g\n", Egrid);
  printf("Etrue   = %g\n", etrue);

  const double deriv_double = f.derive(n);
  Grid grad(nr);
  grad.setZero();
  f.integralgrad(gd, nr, &grad);
  const double deriv_grid = grad[0]/gd.dvolume;
  printf("deriv double = %g\n", deriv_double);
  printf("deriv grid   = %g\n", deriv_grid);

  printf("fractional discrepancy = %g\n", (Edouble - Egrid)/fabs(Edouble));
  if (!(fabs((Edouble - Egrid)/Edouble) < fraccuracy)) {
    printf("FAIL: Discrepancy in the energy is too big!\n");
    retval++;
  }

  printf("fractional error = %g\n", (Edouble - etrue)/fabs(etrue));
  if (!(fabs(Edouble/etrue - 1) < fraccuracy)) {
    printf("FAIL: Error in the energy is too big!\n");
    retval++;
  }

  if (deriv_double == 0) {
    printf("absolute error = %g\n", deriv_double - deriv_grid);
    if (fabs(deriv_double - deriv_grid) > fraccuracy) {
      printf("FAIL: Error in the gradient is too big!\n");
      retval++;
    }
  } else {
    printf("fractional error = %g\n", (deriv_double - deriv_grid)/fabs(deriv_double));
    if (!(fabs((deriv_double - deriv_grid)/deriv_double) < fraccuracy)) {
      printf("FAIL: Error in the gradient is too big!\n");
      retval++;
    }
  }
}

int main(int, char **argv) {
  const double kT = 1e-3;
  const Functional n = EffectivePotentialToDensity(kT);

  {
    double R = 2.2;
    double epsdis = 100;
    double lambda = 1.8;
    const double n = 5e-3;
    const double sphere_volume = 4*M_PI*R*R*R/3;
    const double eta = n*sphere_volume;
    const double c1 = 2.25855 - 1.50349*lambda + 0.249434*lambda*lambda;
    const double c2 = -0.669270 + 1.40049*lambda - 0.827739*lambda*lambda;
    const double c3 = 10.1576 - 15.0427*lambda + 5.30827*lambda*lambda;
    const double eta_effective = c1*eta + c2*eta*eta + c3*eta*eta*eta;
    printf("cn %g %g %g\n", c1, c2, c3);
    printf("eta is %g\n", eta);
    printf("eta_effective is %g\n", eta_effective);
    //const double a1vdw = -4*eta*epsdis*(lambda*lambda*lambda-1);
    const double gHShere = (1-eta_effective/2)/(1-eta_effective)/(1-eta_effective)/(1-eta_effective);
    test_functional("A1", DispersionSAFTa1(R, kT, epsdis, lambda), n,
                    //n*a1vdw*
                    gHShere);

    test_functional("gHScarnahan", gHScarnahan(Identity()*sphere_volume, R), n*eta_effective/eta,
                    gHShere);
    //test_functional("gHS", gHS(Identity()*sphere_volume, R), n, gHShere);
  }

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
