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
#include "LineMinimizer.h"
#include <unistd.h>

#include "utilities.h"

void took(const char *action) {
  static clock_t start = 0;
  clock_t end = clock();
  printf("\n    %s took %g seconds.\n\n", action, (end - double(start))/CLOCKS_PER_SEC);
  start = end;
}

// Here we set up the lattice.
const double zmax = 80;
const double width = 0.0001;
const double cavitysize = 40;

double notinwall(Cartesian r) {
  const double z = r.dot(Cartesian(0,0,1));
  if (fabs(z) > cavitysize/2) {
    return 0;
  } else {
    return 1;
  }
}

int errors = 0;

void test_fast_vs_slow(const char *name, Functional ffast, Functional fslow) {
  
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");


  Lattice lat(Cartesian(width,0,0), Cartesian(0,width,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, 0.1);

  Grid constraint(gd);
  constraint.Set(notinwall);
  //f = constrain(constraint, f);

  Grid potential(gd);

  potential = water_prop.liquid_density*constraint
    + water_prop.vapor_density*VectorXd::Ones(gd.NxNyNz);
  //potential = water_prop.liquid_density*VectorXd::Ones(gd.NxNyNz);

  // Change from density to effective potential...
  potential = -water_prop.kT*potential.cwise().log();
  Grid fastpotential = potential;
  Grid slowpotential = potential;

  Minimizer min = ConjugateGradient(ffast, gd, water_prop.kT, &potential, QuadraticLineMinimizer);

  Minimizer minfast = Precision(1e-18, ConjugateGradient(ffast, gd, water_prop.kT, &fastpotential,
						 QuadraticLineMinimizer));
  Minimizer minslow = Precision(1e-18, ConjugateGradient(fslow, gd, water_prop.kT, &slowpotential,
						 QuadraticLineMinimizer));

  took("Initializing potentials");

  const int num_comparisons = 1;
  for (int i=0;i<num_comparisons && min.improve_energy(true);i++) {
    fflush(stdout);

    double Eslow = fslow.integral(water_prop.kT, potential);
    double Efast = ffast.integral(water_prop.kT, potential);
    double discrepancy = Efast/Eslow - 1;
    printf("Energy discrepancy is %g\n", discrepancy);
    if (fabs(discrepancy) > 1e-12) {
      printf("FAIL: discrepancy in energies is too large!\n");
      errors++;
    }

    Grid gfast(gd), gslow(gd);
    gfast.setZero(); // must zero the gradient before using it!!!
    gslow.setZero(); // must zero the gradient before using it!!!
    fslow.integralgrad(water_prop.kT, potential, &gslow);
    ffast.integralgrad(water_prop.kT, potential, &gfast);
    double gdiscrepancy = (gfast - gslow).cwise().abs().maxCoeff()
      / gfast.cwise().abs().maxCoeff();
    printf("Gradient discrepancy is %g\n", gdiscrepancy);
    if (fabs(gdiscrepancy) > 1e-9) {
      printf("FAIL: discrepancy in gradients is too large!\n");
      errors++;
    }
    // Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    // density.epsNative1d("comparison-water.eps",
		// 	Cartesian(0,0,-zmax/2), Cartesian(0,0,zmax/2),
		// 	water_prop.liquid_density, water_prop.lengthscale, "Y axis: n/n_bulk, x axis: z/R");
    // char *buf = new char[1024];
    // sprintf(buf, "fast-grad-%s.eps", name);
    // gfast.epsNative1d(buf, Cartesian(0,0,-zmax/2), Cartesian(0,0,zmax/2),
    //                   water_prop.liquid_density, water_prop.lengthscale,
    //                   "Y axis: n/n_bulk, x axis: z/R");
    // sprintf(buf, "slow-grad-%s.eps", name);
    // gslow.epsNative1d(buf, Cartesian(0,0,-zmax/2), Cartesian(0,0,zmax/2),
    //                   water_prop.liquid_density, water_prop.lengthscale,
    //                   "Y axis: n/n_bulk, x axis: z/R");
    // delete[] buf;
    // sleep(2);
  }
  took("Verifying energies and gradients match");
  const int numiters = 40;

  Functional myX = Xassociation(water_prop.lengthscale, water_prop.epsilonAB, water_prop.kappaAB,
                                water_prop.epsilon_dispersion, water_prop.lambda_dispersion,
                                water_prop.length_scaling);
  Functional myDelta = DeltaSAFT(water_prop.lengthscale, water_prop.epsilonAB, water_prop.kappaAB,
                                 water_prop.epsilon_dispersion, water_prop.lambda_dispersion,
                                 water_prop.length_scaling);
  printf("\n\nAbout to start improving fast free energy...\n\n");
  for (int i=0;i<numiters && minfast.improve_energy(true);i++) {
    fflush(stdout);
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, fastpotential));
    char *buf = new char[1024];
    sprintf(buf, "fast-%s.eps", name);
    density.epsNative1d(buf, Cartesian(0,0,-zmax/2), Cartesian(0,0,zmax/2),
                        water_prop.liquid_density, water_prop.lengthscale,
                        "Y axis: n/n_bulk, x axis: z/R");

    Grid x(gd, myX(water_prop.kT, gd, density));
    sprintf(buf, "fast-X-%s.eps", name);
    x.epsNative1d(buf, Cartesian(0,0,-zmax/2), Cartesian(0,0,zmax/2),
                  1, water_prop.lengthscale,
                  "y axis: X (unitless), x axis: z/R");

    Grid delta(gd, myDelta(water_prop.kT, gd, density));
    sprintf(buf, "fast-Delta-%s.eps", name);
    delta.epsNative1d(buf, Cartesian(0,0,-zmax/2), Cartesian(0,0,zmax/2),
                      1/water_prop.liquid_density, water_prop.lengthscale,
                      "y axis: delta*nliquid, x axis: z/R");

    Functional n2 = ShellConvolve(water_prop.lengthscale);
    Functional n0 = n2/(4*M_PI*sqr(Functional(water_prop.lengthscale)));
    Grid n0val(gd, n0(water_prop.kT, gd, density));
    sprintf(buf, "fast-n0-%s.eps", name);
    n0val.epsNative1d(buf, Cartesian(0,0,-zmax/2), Cartesian(0,0,zmax/2),
                      water_prop.liquid_density, water_prop.lengthscale,
                      "y axis: n0/nliquid, x axis: z/R");
    delete[] buf;
    // sleep(1);
  }
  took("Minimizing with fast functional");
  printf("\n\nAbout to start improving slow free energy...\n\n");
  for (int i=0;i<numiters && minslow.improve_energy(true);i++) {
    fflush(stdout);
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, slowpotential));
    char *buf = new char[1024];
    sprintf(buf, "slow-%s.eps", name);
    density.epsNative1d(buf, Cartesian(0,0,-zmax/2), Cartesian(0,0,zmax/2),
                        water_prop.liquid_density, water_prop.lengthscale,
                        "Y axis: n/n_bulk, x axis: z/R");
    delete[] buf;
  }
  took("Minimizing with slow functional");

  double fastenergy = minfast.energy()/width/width;
  printf("\nFast energy is %.15g\n", fastenergy);
  double slowenergy = minslow.energy()/width/width;
  printf("Slow energy is %.15g\n", slowenergy);

  printf("Discrepancy in energy is %g\n", fastenergy/slowenergy - 1);
  if (fabs(fastenergy/slowenergy - 1) > 0.1) {
    // FIXME:  We can do better than this!
    printf("FAIL: Discrepancy in energy is too great!\n");
    errors++;
  }

  double Nfast = 0;
  {
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, fastpotential));
    for (int i=0;i<gd.NxNyNz;i++) Nfast += density[i]*gd.dvolume;
  }

  double Nslow = 0;
  {
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, slowpotential));
    for (int i=0;i<gd.NxNyNz;i++) Nslow += density[i]*gd.dvolume;
  }
  printf("Discrepancy in N is %g\n", Nfast/Nslow - 1);
  if (fabs(Nfast/Nslow - 1) > 0.1) {
    // FIXME:  We can do better than this!
    printf("FAIL: Discrepancy in N is too great! (%.15g and %.15g)\n", Nfast, Nslow);
    errors++;
  }

  errors += ffast.run_finite_difference_test("ffast", water_prop.kT, fastpotential);
  errors += fslow.run_finite_difference_test("fslow", water_prop.kT, fastpotential);
}

int main(int, char **argv) {
  Functional f = OfEffectivePotential(SaftFluid(water_prop.lengthscale,
                                                water_prop.epsilonAB, water_prop.kappaAB,
                                                water_prop.epsilon_dispersion,
                                                water_prop.lambda_dispersion, water_prop.length_scaling, 0));
  double mu = find_chemical_potential(f, water_prop.kT,
                                      1.01*water_prop.liquid_density);
  Functional fslow = OfEffectivePotential(SaftFluidSlow(water_prop.lengthscale,
                                                        water_prop.epsilonAB, water_prop.kappaAB,
                                                        water_prop.epsilon_dispersion,
                                                        water_prop.lambda_dispersion, water_prop.length_scaling, mu));
  Functional ffast = OfEffectivePotential(SaftFluid(water_prop.lengthscale,
                                                    water_prop.epsilonAB, water_prop.kappaAB,
                                                    water_prop.epsilon_dispersion,
                                                    water_prop.lambda_dispersion, water_prop.length_scaling, mu));
  took("Creating functionals");

  test_fast_vs_slow("Saft", ffast, fslow);

  Functional n = EffectivePotentialToDensity();
  f = HardSpheresWBnotensor(water_prop.lengthscale)(n);
  mu = find_chemical_potential(f, water_prop.kT, 1.3*water_prop.liquid_density);
  fslow = OfEffectivePotential(HardSpheresWBnotensor(water_prop.lengthscale)
                               + ChemicalPotential(mu));
  ffast = OfEffectivePotential(HardSpheresNoTensor(water_prop.lengthscale)
                               + ChemicalPotential(mu));

  test_fast_vs_slow("HardSpheresNoTensor", ffast, fslow);

  if (errors == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], errors);
  return errors;
}
