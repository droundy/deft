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
#include <time.h>
#include "Functionals.h"
#include "utilities.h"

int retval = 0;
double kT = water_prop.kT;

void test_eos(const char *name, Functional f, double ntrue, double ptrue, double fraccuracy=1e-6) {
  clock_t start = clock();
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  printf("Expect Veff of %g\n", -kT*log(ntrue));
  printf("Looking for density between %g and %g...\n", ntrue*3.123e-7, ntrue*12345);
  double nfound = find_density(f, kT, ntrue*3.123e-7, ntrue*12345);
  printf("Found density of %.15g (versus %g) in %g seconds.\n", nfound, ntrue,
         (clock() - double(start))/CLOCKS_PER_SEC);
  double nerr = nfound/ntrue - 1;
  if (fabs(nerr) > fraccuracy) {
    printf("FAIL: nerr too big: %g\n", nerr);
    retval++;
  }

  double pfound = pressure(f, kT, nfound);
  printf("Found pressure of %.15g (versus %g) in %g seconds.\n", pfound, ptrue,
         (clock() - double(start))/CLOCKS_PER_SEC);
  double perr = pfound/ptrue - 1;
  if (fabs(perr) > fraccuracy) {
    printf("FAIL: pnerr too big: %g\n", perr);
    retval++;
  }
}

void test_pressure(const char *name, Functional f, double n, double ptrue, double fraccuracy=1e-6) {
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  double pfound = pressure(f, kT, n);
  printf("Found pressure of %.15g (versus %g).\n", pfound, ptrue);
  double perr = pfound/ptrue - 1;
  if (fabs(perr) > fraccuracy) {
    printf("FAIL: pnerr too big: %g\n", perr);
    retval++;
  }
}

int main(int, char **argv) {
  Functional n = EffectivePotentialToDensity(kT);
  double Veff = -kT*log(water_prop.liquid_density);

  {
    double ngas = 2e-5;
    double mu = -kT*log(ngas);
    test_eos("ideal gas", IdealGasOfVeff(kT) + ChemicalPotential(mu)(n), ngas, ngas*kT);
  }

  test_eos("quadratic", 0.5*sqr(n) - n, 1.0, 0.5, 2e-6);
  test_pressure("quadratic(2)", 0.5*sqr(n) - n, 2, 2);
  test_pressure("quadratic(3)", 0.5*sqr(n) - n, 3, 4.5);

  {
    FILE *o = fopen("ideal-gas.dat", "w");
    equation_of_state(o, IdealGasOfVeff(kT), kT, 1e-7, 1e-3);
    fclose(o);
  }

  {
    FILE *o = fopen("saft-fluid.dat", "w");
    Functional f = SaftFluidSlow(water_prop.lengthscale, kT,
                                 water_prop.epsilonAB, water_prop.kappaAB,
                                 water_prop.epsilon_dispersion, water_prop.lambda_dispersion, 0);
    double mu = f.derive(Veff)*kT/water_prop.liquid_density; // convert from derivative w.r.t. V
    equation_of_state(o, f + ChemicalPotential(mu)(n), kT, 1e-7, 7e-3);
    fclose(o);

    o = fopen("saft-fluid-other.dat", "w");
    //other_equation_of_state(o, f + ChemicalPotential(mu)(n), kT, 1e-7, 7e-3);
    fclose(o);

    Functional X = Xassociation(water_prop.lengthscale, kT,
                                water_prop.epsilonAB, water_prop.kappaAB);
    printf("X is %g\n", X(water_prop.liquid_density));
    o = fopen("association.dat", "w");
    equation_of_state(o, AssociationSAFT(water_prop.lengthscale, kT,
                                         water_prop.epsilonAB, water_prop.kappaAB)(n),
                      kT, 1e-7, 0.0095);
    fclose(o);
    o = fopen("dispersion.dat", "w");
    equation_of_state(o, DispersionSAFT(water_prop.lengthscale, kT,
                                        water_prop.epsilon_dispersion,
                                        water_prop.lambda_dispersion)(n),
                      kT, 1e-7, 0.0095);
    fclose(o);
  }

  {
    FILE *o = fopen("hard-sphere-fluid.dat", "w");
    Functional f = HardSpheresWBnotensor(water_prop.lengthscale, kT)(n)
      + IdealGasOfVeff(kT);
    double mu = f.derive(Veff)*kT/water_prop.liquid_density; // convert from derivative w.r.t. V
    equation_of_state(o, f + ChemicalPotential(mu)(n), kT, 1e-7, 7e-3);
    fclose(o);
  }

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
