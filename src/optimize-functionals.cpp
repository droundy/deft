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

#include "Functionals.h"
#include "utilities.h"
#include <stdlib.h>

const double R = 2.7;
const double mu = 0.01;

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("Please provide a filename\n");
    exit(1);
  }
  const char *R_arg[] = { "R", 0 };
  const char *R_mu_arg[] = { "R", "mu", 0 };
  if (strcmp(argv[1], "src/HardSpheresFast.cpp") == 0)
    HardSpheres(R).create_source("src/HardSpheresFast.cpp", "HardSpheresFast", R_arg);
  if (strcmp(argv[1], "src/HardSpheresRFFast.cpp") == 0)
    HardSpheresRF(R).create_source("src/HardSpheresRFFast.cpp", "HardSpheresRFFast", R_arg);
  if (strcmp(argv[1], "src/HardSpheresTarazonaFast.cpp") == 0)
    HardSpheresTarazona(R).create_source("src/HardSpheresTarazonaFast.cpp", "HardSpheresTarazonaFast", R_arg);
  if (strcmp(argv[1], "src/HardSpheresNoTensorFast.cpp") == 0)
    HardSpheresWBnotensor(R).create_source(argv[1], "HardSpheresNoTensor", R_arg);
  if (strcmp(argv[1], "src/HardSpheresWBm2Fast.cpp") == 0)
    HardSpheresWBm2slow(R).create_source(argv[1], "HardSpheresWBm2", R_arg);
  if (strcmp(argv[1], "src/HardSpheresWBFast.cpp") == 0)
    HardSpheresWB(R).create_source(argv[1], "HardSpheresWBFast", R_arg);
  if (strcmp(argv[1], "src/HardSphereGasRFFast.cpp") == 0) {
    Functional f = HardSpheresRF(R) + IdealGas() + ChemicalPotential(mu);
    f.create_source(argv[1], "HardSphereGasRF", R_mu_arg);
  }
  if (strcmp(argv[1], "src/HardSphereGasFast.cpp") == 0) {
    Functional f = HardSpheres(R) + IdealGas() + ChemicalPotential(mu);
    f.create_source(argv[1], "HardSphereGas", R_mu_arg);
  }
  const char *saft_args[] = {
    "R",
    "epsilonAB", "kappaAB",
    "epsilon_dispersion", "lambda_dispersion", "length_scaling", "mu",
    0
  };
  if (strcmp(argv[1], "src/SaftFluidFast.cpp") == 0)
    SaftFluidSlow(R, 0.01, 0.01, 0.01, 0.01, 0.01,
                  mu).create_source(argv[1], "SaftFluid", saft_args);
  const char *dispersion_args[] = { "R", "epsilon_dispersion", "lambda_dispersion",
                                    "length_scaling", 0 };
  if (strcmp(argv[1], "src/DispersionFast.cpp") == 0) {
    DispersionSAFT(R, 0.01, 0.01, 0.01).create_source(argv[1], "Dispersion", dispersion_args);
  }
  const char *association_args[] = {
    "R",
    "epsilonAB", "kappaAB",
    "epsilon_dispersion", "lambda_dispersion", "length_scaling",
    0
  };
  if (strcmp(argv[1], "src/AssociationFast.cpp") == 0)
    AssociationSAFT(R, 0.01, 0.01, 0.01,
                    0.01, 0.01).create_source(argv[1], "Association", association_args);
}
