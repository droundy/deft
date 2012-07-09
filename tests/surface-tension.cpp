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
#include "handymath.h"

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  printf("%s took %g seconds\n", name, (t-last_time)/double(CLOCKS_PER_SEC));
  last_time = t;
}

const double nv = 1e-7, nl = 1e-3;
Functional n0 = Gaussian(1);
// The following is the simplest interesting functional
Functional f = /*IdealGasOfVeff() + */ 1e6*sqr(n0 - Functional(water_prop.liquid_density))*sqr(n0 - Functional(water_prop.vapor_density));

int main(int, char **argv) {
  int retval = 0;
  const double kB = 3.16681539628059e-6; // Boltzmann's constant in Hartree/Kelvin
  LiquidProperties prop = water_prop;
  prop.lengthscale = 2;
  double criten = 0;
  for (double n = prop.vapor_density; n < prop.liquid_density; n *= 1.1) {
    if (f(prop.kT, n) > criten) {
      criten = f(prop.kT, n);
      prop.critical_density = n;
    }
    // printf("f(%g) = %g\n", n, f(prop.kT, n));
  }
  printf("Critical density: %g\n", prop.critical_density);
  f = OfEffectivePotential(f);
  saturated_liquid_properties(f, &prop);
  prop.critical_density = 1;
  took("Finding bulk densities");
  if (fabs(prop.liquid_density-water_prop.liquid_density) > 1e-9*water_prop.liquid_density) {
    printf("FAIL: bad liquid density:  %g (vs. %g)\n", prop.liquid_density, water_prop.liquid_density);
    retval++;
  }
  if (fabs(prop.vapor_density-water_prop.vapor_density) > 1e-9*water_prop.vapor_density) {
    printf("FAIL: bad vapor density:  %g (vs. %g)\n", prop.vapor_density, water_prop.vapor_density);
    retval++;
  }
  // Here we set up an unused lattice.
  Lattice lat(Cartesian(0.2,0,0), Cartesian(0,0.2,0), Cartesian(0,0,20));
  GridDescription gd(lat, 1, 1, 200);
  Grid foo(gd);
  Minimizer min = MaxIter(500, Precision(1e-14, ConjugateGradient(f, gd, prop.kT, &foo, QuadraticLineMinimizer)));
  const bool amverbose = false;
  const double st = surface_tension(min, f, prop, amverbose);
  took("Finding surface tension");
  const double true_st = 4.08612e-5;
  printf("surface tension is %.15g (vs %.15g) (%g fractionally)\n", st, true_st, fabs(st-true_st)/true_st);
  if (fabs(st - true_st) > 1e-3*true_st) {
    printf("FAIL: bad surface tension: %.15g (vs %.15g) (%g fractionally)\n", st, true_st, fabs(st-true_st)/true_st);
    retval++;
  }
  if (retval > 0) {
    printf("%s failed %d tests!\n", argv[0], retval);
  } else {
    printf("%s passed!\n", argv[0]);
  }
  return retval;
}
