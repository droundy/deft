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
#include "new/WaterSaftFast.h"
#include "new/HomogeneousWaterSaftFast.h"
#include "new/WaterSaftByHandFast.h"
#include "new/HomogeneousWaterSaftByHandFast.h"
#include "equation-of-state.h"

int check_functional_value(const char *name,
                            const NewFunctional &f,
                            double energy,
                            double fraccuracy = 1e-15) {
  int errors = 0;
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  double fv = f.energy();
  print_double("Energy of Vector:  ", fv);
  printf("\n");
  f.printme("  ");

  if (!(fabs(fv/energy - 1) < fraccuracy)) {
    printf("fv = %g\n", fv);
    printf("expected = %g\n", energy);
    printf("FAIL: Error in f(n) is %g\n", fv/energy - 1);
    errors++;
  }
  return errors;
}

int main(int, char **argv) {
  int retval = 0;

  const int Nx = 100;
  const double a = 10.0, kT = new_water_prop.kT, nval = 1.1*new_water_prop.liquid_density;
  const double energy_density = -9.350583155126999e-05; // from older code
  printf("about to create input\n");
  WaterSaft ws(Nx, Nx, Nx);

  // double lengthscale, liquid_density, critical_density, vapor_density, kT;
  // double epsilon_dispersion, lambda_dispersion, length_scaling, epsilonAB, kappaAB;

  ws.R() = new_water_prop.lengthscale;
  ws.epsilon_association() = new_water_prop.epsilonAB;
  ws.epsilon_dispersion() = new_water_prop.epsilon_dispersion;
  ws.kappa_association() = new_water_prop.kappaAB;
  ws.lambda_dispersion() = new_water_prop.lambda_dispersion;
  ws.length_scaling() = new_water_prop.length_scaling;
  ws.a1() = a;
  ws.a2() = a;
  ws.a3() = a;
  ws.kT() = kT;
  ws.n() = nval;
  retval += check_functional_value("WaterSaft", ws, energy_density*a*a*a, 1e-11);

  WaterSaftByHand wsbh(Nx, Nx, Nx);

  // double lengthscale, liquid_density, critical_density, vapor_density, kT;
  // double epsilon_dispersion, lambda_dispersion, length_scaling, epsilonAB, kappaAB;

  wsbh.R() = new_water_prop.lengthscale;
  wsbh.epsilon_association() = new_water_prop.epsilonAB;
  wsbh.epsilon_dispersion() = new_water_prop.epsilon_dispersion;
  wsbh.kappa_association() = new_water_prop.kappaAB;
  wsbh.lambda_dispersion() = new_water_prop.lambda_dispersion;
  wsbh.length_scaling() = new_water_prop.length_scaling;
  wsbh.a1() = a;
  wsbh.a2() = a;
  wsbh.a3() = a;
  wsbh.kT() = kT;
  wsbh.n() = nval;
  retval += check_functional_value("WaterSaftByHand", wsbh, energy_density*a*a*a, 1e-11);

  HomogeneousWaterSaft hws;
  hws.R() = new_water_prop.lengthscale;
  hws.epsilon_association() = new_water_prop.epsilonAB;
  hws.epsilon_dispersion() = new_water_prop.epsilon_dispersion;
  hws.kappa_association() = new_water_prop.kappaAB;
  hws.lambda_dispersion() = new_water_prop.lambda_dispersion;
  hws.kT() = kT;
  hws.n() = nval;
  retval += check_functional_value("HomogeneousWaterSaft", hws, energy_density);

  HomogeneousWaterSaftByHand hwsbh;
  hwsbh.R() = new_water_prop.lengthscale;
  hwsbh.epsilon_association() = new_water_prop.epsilonAB;
  hwsbh.epsilon_dispersion() = new_water_prop.epsilon_dispersion;
  hwsbh.kappa_association() = new_water_prop.kappaAB;
  hwsbh.lambda_dispersion() = new_water_prop.lambda_dispersion;
  hwsbh.kT() = kT;
  hwsbh.n() = nval;
  retval += check_functional_value("HomogeneousWaterSaftByHand", hwsbh, energy_density);

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
