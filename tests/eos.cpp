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
#include "equation-of-state.h"

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

  const double nmin = 1e-11, nmax = 0.007;

  {
    double ngas = 2e-5;
    double mu = find_chemical_potential(IdealGasOfVeff(kT), kT, ngas);
    test_eos("ideal gas", IdealGasOfVeff(kT) + ChemicalPotential(mu)(n), ngas, ngas*kT);
  }

  test_eos("quadratic", 0.5*sqr(n) - n, 1.0, 0.5, 2e-6);
  test_pressure("quadratic(2)", 0.5*sqr(n) - n, 2, 2);
  test_pressure("quadratic(3)", 0.5*sqr(n) - n, 3, 4.5);

  {
    //FILE *o = fopen("ideal-gas.dat", "w");
    //equation_of_state(o, IdealGasOfVeff(kT), kT, nmin, nmax);
    //fclose(o);
  }

  {
    FILE *o = fopen("dispersion.dat", "w");
    //equation_of_state(o, DispersionSAFT(water_prop.lengthscale, kT,
    //                                    water_prop.epsilon_dispersion,
    //                                    water_prop.lambda_dispersion)(n),
    //                  kT, nmin, nmax);
    fclose(o);
    printf("Got dispersion!\n");

    Functional f = SaftFluidSlow(water_prop.lengthscale, kT,
                                 water_prop.epsilonAB, water_prop.kappaAB,
                                 water_prop.epsilon_dispersion,
                                 water_prop.lambda_dispersion, 0);

    {
      double nv = coexisting_vapor_density(f, water_prop.kT, water_prop.liquid_density);
      printf("predicted vapor density: %g\n", nv);
      printf("actual vapor density:    %g\n", water_prop.vapor_density);
    }

    if (0) {
      o = fopen("saft-fluid.dat", "w");
      double mu = f.derive(Veff)*kT/water_prop.liquid_density; // convert from derivative w.r.t. V
      equation_of_state(o, f + ChemicalPotential(mu)(n), kT, nmin, nmax);
      fclose(o);
    }

    {
      const double nl = saturated_liquid(f, water_prop.kT);
      printf("saturated water density is %g\n", nl);
      if (fabs(nl/water_prop.liquid_density - 1) > 0.1) {
        printf("FAIL: error in saturated water density is too big! %g\n",
               nl/water_prop.liquid_density - 1);
        retval++;
      }

      double nv = coexisting_vapor_density(f, water_prop.kT, nl);
      printf("predicted saturated vapor density: %g\n", nv);
      printf("actual vapor density:    %g\n", water_prop.vapor_density);
      //double mu = f.derive(-kT*log(nl))*kT/nl; // convert from derivative w.r.t. V
      //o = fopen("saft-fluid-saturated.dat", "w");
      //equation_of_state(o, f + ChemicalPotential(mu)(n), kT, nmin, 1.1*nl);
      //fclose(o);

      double pv = pressure(f, water_prop.kT, nv);
      printf("vapor pressure is %g\n", pv);
      if (fabs(pv/water_prop.kT/nv - 1) > 1e-3) {
        printf("FAIL: error in vapor pressure, steam isn't ideal gas? %g\n",
               pv/water_prop.kT/nv - 1);
        retval++;
      }
    }

    {
      o = fopen("room-temperature.dat", "w");
      Functional f = SaftFluidSlow(water_prop.lengthscale, water_prop.kT,
                                   water_prop.epsilonAB, water_prop.kappaAB,
                                   water_prop.epsilon_dispersion,
                                   water_prop.lambda_dispersion, 0);
      double mufoo = find_chemical_potential(f, water_prop.kT,
                                             water_prop.liquid_density);
      f = SaftFluidSlow(water_prop.lengthscale, kT,
                        water_prop.epsilonAB, water_prop.kappaAB,
                        water_prop.epsilon_dispersion,
                        water_prop.lambda_dispersion, mufoo);
      printf("moofoo is %g\n", mufoo);
      double nl = saturated_liquid(f, water_prop.kT);
      double nv = coexisting_vapor_density(f, water_prop.kT, nl);
      //double p = pressure(f, water_prop.kT, nl);
      double mu = find_chemical_potential(f, water_prop.kT, nl);
      for (double dens=0.1*nv; dens<=1.2*nl; dens *= 1.01) {
        double V = -water_prop.kT*log(dens);
        double Vl = -water_prop.kT*log(nl);
        fprintf(o, "%g\t%g\t%g\n",
                dens, f(V), f(Vl) - (dens-nl)*mu);
      }
      fclose(o);
      printf("Finished plotting room-temperature.dat...\n");
    }

    const double n_1atm = pressure_to_density(f, water_prop.kT, atmospheric_pressure);
    printf("density at 1 atmosphere is %g\n", n_1atm);
    if (fabs(n_1atm/water_prop.liquid_density - 1) > 0.1) {
      printf("FAIL? error in water density is too big! %g\n",
             n_1atm/water_prop.liquid_density - 1);
      //retval++;
    }

    /*
    o = fopen("saft-fluid-other.dat", "w");
    //other_equation_of_state(o, f + ChemicalPotential(mu)(n), kT, 1e-7, 7e-3);
    fclose(o);

    Functional X = Xassociation(water_prop.lengthscale, kT,
                                water_prop.epsilonAB, water_prop.kappaAB);
    printf("X is %g\n", X(water_prop.liquid_density));
    */
    o = fopen("association.dat", "w");
    equation_of_state(o, AssociationSAFT(water_prop.lengthscale, kT,
                                         water_prop.epsilonAB, water_prop.kappaAB,
                                         water_prop.epsilon_dispersion,
                                         water_prop.lambda_dispersion)(n),
                      kT, nmin, nmax);
    fclose(o);
  }

  {
    FILE *o = fopen("hard-sphere-fluid.dat", "w");
    Functional f = HardSpheresWBnotensor(water_prop.lengthscale, kT)(n)
      + IdealGasOfVeff(kT);
    double mu = f.derive(Veff)*kT/water_prop.liquid_density; // convert from derivative w.r.t. V
    equation_of_state(o, f + ChemicalPotential(mu)(n), kT, nmin, nmax);
    fclose(o);
  }

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
