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
#include "equation-of-state.h"

int retval = 0;

const double kT = 1e-3; // room temperature in Hartree, approximately

void test_energy(const char *name, Functional f,
                 double true_energy, double fraccuracy = 1e-15) {
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  Lattice lat(Cartesian(0,.5,.5), Cartesian(.5,0,.5), Cartesian(.5,.5,0));
  int resolution = 20;
  GridDescription gd(lat, resolution, resolution, resolution);
  Grid density(gd);
  density = 1e-3*(-500*r2(gd)).cwise().exp()
    + 1e-7*VectorXd::Ones(gd.NxNyNz);
  Grid eff_potential(gd, -kT*density.cwise().log());

  retval += f.run_finite_difference_test(name, eff_potential);

  double e = f.integral(eff_potential);
  printf("Energy = %.16g\n", e);
  printf("Fractional error = %g\n", (e - true_energy)/fabs(true_energy));
  if (!(fabs((e - true_energy)/true_energy) < fraccuracy)) {
    printf("FAIL: Error in the energy is too big!\n");
    retval++;
  }
}

int main(int, char **argv) {
  Functional n = EffectivePotentialToDensity(kT);
  test_energy("association",
              AssociationSAFT(water_prop.lengthscale, kT,
                              water_prop.epsilonAB, water_prop.kappaAB,
                              water_prop.epsilon_dispersion,
                              water_prop.lambda_dispersion)(n),
              -4.648379504003813e-12);
  const double dispersion_energy = -2.255139236942726e-12;
  test_energy("dispersion",
              DispersionSAFT(water_prop.lengthscale, kT,
                             water_prop.epsilon_dispersion,
                             water_prop.lambda_dispersion)(n),
              dispersion_energy);
  {
    // The following should work, once I properly split up the
    // dispersion energy...
    Functional R(water_prop.lengthscale);
    R.set_name("R");
    Functional n2 = ShellConvolve(water_prop.lengthscale);
    Functional n0 = n2/(4*M_PI*sqr(R));
    Functional a1 = DispersionSAFTa1(water_prop.lengthscale,
                                     water_prop.epsilon_dispersion,
                                     water_prop.lambda_dispersion);
    Functional a2 = DispersionSAFTa2(water_prop.lengthscale,
                                     water_prop.epsilon_dispersion,
                                     water_prop.lambda_dispersion);
    // The following is bogus because I no longer scale by n0...
    //test_energy("dispersion by parts", (n0*(a1 + a2/kT))(n), dispersion_energy);
  }
  test_energy("SAFT slow",
              SaftFluidSlow(water_prop.lengthscale, kT,
                            water_prop.epsilonAB, water_prop.kappaAB,
                            water_prop.epsilon_dispersion,
                            water_prop.lambda_dispersion, 0),
             -8.140476801563833e-09);
  //test_energy("SAFT",
  //            SaftFluid(water_prop.lengthscale, kT,
  //                      water_prop.epsilonAB, water_prop.kappaAB,
  //                      water_prop.epsilon_dispersion,
  //                      water_prop.lambda_dispersion, 0),
  //            1.0);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
