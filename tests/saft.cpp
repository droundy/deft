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

void test_energy(const char *name, Functional f, double kT,
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

  retval += f.run_finite_difference_test(name, kT, density);

  double e = f.integral(kT, density);
  printf("Energy = %.16g\n", e);
  printf("Fractional error = %g\n", (e - true_energy)/fabs(true_energy));
  if (e < true_energy) {
    printf("FAIL: the energy is too low! (it is %.16g)\n", e);
    retval++;
  }
  if (!(fabs((e - true_energy)/true_energy) < fraccuracy)) {
    printf("FAIL: Error in the energy is too big!\n");
    retval++;
  }
}

int main(int, char **argv) {
  const double kT = 1e-3; // room temperature in Hartree, approximately
  test_energy("association",
              AssociationSAFT(water_prop.lengthscale,
                              water_prop.epsilonAB, water_prop.kappaAB,
                              water_prop.epsilon_dispersion,
                              water_prop.lambda_dispersion, water_prop.length_scaling),
              kT, -4.663737042961209e-12);
  const double dispersion_energy = -2.25551876171605e-12;
  test_energy("dispersion",
              DispersionSAFT(water_prop.lengthscale,
                             water_prop.epsilon_dispersion,
                             water_prop.lambda_dispersion, water_prop.length_scaling),
              kT, dispersion_energy);
  test_energy("SAFT slow",
              SaftFluidSlow(water_prop.lengthscale,
                            water_prop.epsilonAB, water_prop.kappaAB,
                            water_prop.epsilon_dispersion,
                            water_prop.lambda_dispersion, water_prop.length_scaling, 0),
              kT, -8.140492183697745e-09);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
