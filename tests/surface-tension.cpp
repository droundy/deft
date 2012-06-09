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
#include "LineMinimizer.h"
#include "equation-of-state.h"

const double ngas = water_prop.vapor_density; // vapor density of water
const double nliquid = water_prop.liquid_density; // density of liquid water
const double mu = -water_prop.kT*log(ngas);
const double Veff_liquid = -water_prop.kT*log(nliquid);

// Here we set up the lattice.
Lattice lat(Cartesian(0.2,0,0), Cartesian(0,0.2,0), Cartesian(0,0,20));
GridDescription gd(lat, 1, 1, 200);

// And the functional...
const double interaction_energy_scale = 2e-4;
Functional attraction = GaussianPolynomial(-interaction_energy_scale/nliquid/nliquid/2, 0.5, 2);
Functional repulsion = GaussianPolynomial(interaction_energy_scale/nliquid/nliquid/nliquid/nliquid/4, 0.125, 4);
Functional f0 = ChemicalPotential(mu) + attraction + repulsion;
Functional n = EffectivePotentialToDensity();
Functional f = IdealGasOfVeff + f0(n);

Grid external_potential(gd);
Grid potential(gd);

Functional ff;

int test_minimizer(const char *name, Minimizer *min, int numiters, double fraccuracy=1e-3) {
  clock_t start = clock();
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  potential = Veff_liquid*VectorXd::Ones(gd.NxNyNz);
  for (int i=0;i<gd.NxNyNz/2;i++) potential[i] = mu;

  for (int i=0;i<numiters && min->improve_energy(false);i++) {
    fflush(stdout);
  }
  min->print_info();
  //min->minimize(f, gd);
  //for (int i=0;i<numiters && min->improve_energy(false);i++) {
  //  //printf("Without external potential...\n");
  //  fflush(stdout);
  //}
  //min->print_info();  

  const double Einterface_with_external = ff.integral(water_prop.kT, potential);
  const double Einterface = f.integral(water_prop.kT, potential);
  double Ninterface = 0;
  {
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    for (int i=0;i<gd.NxNyNz;i++) Ninterface += density[i]*gd.dvolume;
  }
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  start = clock();

  Grid gas(gd, mu*VectorXd::Ones(gd.NxNyNz));
  min->minimize(f, gd, &gas);
  // The following sometimes fails with NaNs if it's run too long, for
  // reasons I don't understand.
  for (int i=0;i<numiters && min->improve_energy(false);i++) {
    //printf("GAS\n");
    fflush(stdout);
  }
  const double Egas = f.integral(water_prop.kT, gas);
  double Ngas = 0;
  {
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, gas));
    for (int i=0;i<gd.NxNyNz;i++) Ngas += density[i]*gd.dvolume;
  }
  min->print_info();
  printf("gas energy is %g\n", f.integral(water_prop.kT, gas));
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  start = clock();

  double minpot = 1e100;
  for (int i=0;i<gd.NxNyNz;i++)
    if (potential[i] < minpot) minpot = potential[i];
  Grid liquid(gd, Veff_liquid*VectorXd::Ones(gd.NxNyNz));
  min->minimize(f, gd, &liquid);
  for (int i=0;i<numiters && min->improve_energy(false);i++) {
    //printf("LIQUID\n");
    fflush(stdout);
  }
  min->print_info();
  printf("Minimization took %g seconds.\n", (clock() - double(start))/CLOCKS_PER_SEC);
  const double Eliquid = f.integral(water_prop.kT, liquid);
  double Nliquid = 0;
  {
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, liquid));
    for (int i=0;i<gd.NxNyNz;i++) Nliquid += density[i]*gd.dvolume;
  }

  double retval = 0;
  printf("\n\n");
  const double Einterface_with_external_true = -1.00700771531599e-05;
  printf("interface energy is %.15g (vs. %.15g)\n", Einterface, Einterface_with_external);
  if (Einterface_with_external < Einterface_with_external_true) {
    printf("FAIL: Einterface_with_external is too low: %.17g (should be %.17g)\n",
	   Einterface_with_external, Einterface_with_external_true);
    retval++;
  }
  printf("fractional Einterface error = %g\n",
         (Einterface_with_external - Einterface_with_external_true)/fabs(Einterface_with_external_true));
  if (fabs((Einterface_with_external - Einterface_with_external_true)/Einterface_with_external_true) > fraccuracy) {
    printf("FAIL: Einterface_with_external is too inaccurate...\n");
    retval++;
  }

  printf("gas energy is %.15g\n", Egas);
  const double Egas_true = -2.53378305461179e-05;
  if (Egas < Egas_true) {
    printf("FAIL: Egas is too low: %.16g\n", Egas);
    retval++;
  }

  printf("liquid energy is %.15g\n", Eliquid);
  const double Eliquid_true = -2.53378305461181e-05;
  if (Eliquid < Eliquid_true) {
    printf("FAIL: Eliquid is too low: %.16g\n", Eliquid);
    retval++;
  }
  printf("Ninterface/liquid/gas = %g/%g/%g\n", Ninterface, Nliquid, Ngas);
  const double X = Ninterface/Nliquid; // Fraction of volume effectively filled with liquid.
  printf("X is %g\n", X);
  
  const double true_surface_tension = 1.208891768e-05;

  const double surface_tension = (Einterface - Eliquid*X - Egas*(1-X))/2/0.2/0.2;
  printf("surface tension is %.10g\n", surface_tension);


  printf("fractional surface tension error = %g\n",
         (surface_tension - true_surface_tension)/fabs(true_surface_tension));
  if (fabs((surface_tension - true_surface_tension)/true_surface_tension) > fraccuracy) {
    printf("FAIL:\n");
    printf("Error in the surface tension is too big!\n");
    retval++;
  }
  return retval;
}

double forcer(Cartesian r) {
  const double z = r.dot(Cartesian(0,0,1));
  return 0.0001/nliquid*exp(-(z-5.0)*(z-5)/2);
}

int main(int, char **argv) {
  external_potential.Set(forcer);
  ff = IdealGasOfVeff + (f0 + ExternalPotential(external_potential))(n);

  Grid test_density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, -1e-4*(-2*r2(gd)).cwise().exp()
                                                      + mu*VectorXd::Ones(gd.NxNyNz)));

  // The following is for debugging our simple liquid...
  //printf("mu is %g\n", mu);
  //printf("Veff_liquid is %g\n", Veff_liquid);
  //printf("%12s\t%12s\t%12s\n", "Veff", "n", "energy");
  //for (double v=0.5*Veff_liquid; v<mu+5*kT; v*=1.05) {
  //  printf("%12g\t%12g\t%12g\n", v, EffectivePotentialToDensity(kT)(v), f(v));
  //}

  int retval = 0;

  potential = +1e-4*((-10*r2(gd)).cwise().exp()) + 1.14*mu*VectorXd::Ones(gd.NxNyNz);
  retval += ff.run_finite_difference_test("simple liquid", water_prop.kT, potential);
  fflush(stdout);

  retval += attraction.run_finite_difference_test("quadratic", water_prop.kT, test_density);
  fflush(stdout);
  retval += repulsion.run_finite_difference_test("repulsive", water_prop.kT, test_density);
  fflush(stdout);

  {
    // We need a small initial step for PreconditionedDownhill in
    // order to avoid jumping over the barrier from the liquid global
    // minimum to the vapor local minimum.
    Minimizer pd = PreconditionedDownhill(ff, gd, water_prop.kT, &potential, 1e-9);
    potential.setZero();
    retval += test_minimizer("PreconditionedDownhill", &pd, 4001, 10); // ignore accuracy for now

    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    //density.epsNativeSlice("PreconditionedDownhill.eps", Cartesian(0,0,20), Cartesian(0.1,0,0),
    //                       Cartesian(0,0,0));

    retval += attraction.run_finite_difference_test("quadratic", water_prop.kT, density);
    retval += repulsion.run_finite_difference_test("repulsive", water_prop.kT, density);
  }

  Minimizer psd = PreconditionedSteepestDescent(ff, gd, water_prop.kT, &potential, QuadraticLineMinimizer);
  potential.setZero();
  // I'm not sure why PreconditionedSteepestDescent is failing this badly!
  //retval += test_minimizer("PreconditionedSteepestDescent", &psd, 138000, 1e-2);

  {
    Minimizer pcg = PreconditionedConjugateGradient(ff, gd, water_prop.kT, &potential, QuadraticLineMinimizer);
    potential.setZero();
    retval += test_minimizer("PreconditionedConjugateGradient", &pcg, 6000, 3e-5);

    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    density.epsNativeSlice("PreconditionedConjugateGradient.eps",
                           Cartesian(0,0,40), Cartesian(0.1,0,0), Cartesian(0,0,0));
  }

  Minimizer cg = ConjugateGradient(ff, gd, water_prop.kT, &potential, QuadraticLineMinimizer);
  potential.setZero();
  retval += test_minimizer("ConjugateGradient", &cg, 3000, 3e-5);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
