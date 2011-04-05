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
  if (strcmp(argv[1], "src/IdealGasFast.cpp") == 0)
    IdealGasOfVeff.create_source("src/IdealGasFast.cpp", "IdealGasFast");
  if (strcmp(argv[1], "src/HardSpheresFast.cpp") == 0)
    HardSpheres(R).create_source("src/HardSpheresFast.cpp", "HardSpheresFast", "R");
  if (strcmp(argv[1], "src/HardSpheresRFFast.cpp") == 0)
    HardSpheresRF(R).create_source("src/HardSpheresRFFast.cpp", "HardSpheresRFFast", "R");
  if (strcmp(argv[1], "src/HardSpheresTarazonaFast.cpp") == 0)
    HardSpheresTarazona(R).create_source("src/HardSpheresTarazonaFast.cpp", "HardSpheresTarazonaFast", "R");
  if (strcmp(argv[1], "src/HardSpheresNoTensorFast.cpp") == 0)
    HardSpheresWBnotensor(R).create_source(argv[1], "HardSpheresNoTensor", "R");
  if (strcmp(argv[1], "src/HardSphereGasRFFast.cpp") == 0) {
    Functional n = EffectivePotentialToDensity();
    Functional f = HardSpheresRF(R)(n) + IdealGasOfVeff + ChemicalPotential(mu)(n);
    f.create_source(argv[1], "HardSphereGasRF", "R", "mu");
  }
  if (strcmp(argv[1], "src/HardSphereGasFast.cpp") == 0) {
    Functional n = EffectivePotentialToDensity();
    Functional f = HardSpheres(R)(n) + IdealGasOfVeff + ChemicalPotential(mu)(n);
    f.create_source(argv[1], "HardSphereGas", "R", "mu");
  }
  if (strcmp(argv[1], "src/SaftFluidFast.cpp") == 0) {
    SaftFluidSlow(R, 0.01, 0.01, 0.01, 0.01, mu).create_source(argv[1], "SaftFluid", "R",
                                                               "epsilonAB", "kappaAB",
                                                               "epsilon_dispersion", "lambda_dispersion",
                                                               "mu");
  }
  if (strcmp(argv[1], "src/SaftExcessEnergyFast.cpp") == 0) {
    SaftExcessEnergySlow(R, 0.01, 0.01, 0.01, 0.01,
                         mu).create_source(argv[1], "SaftExcessEnergy", "R",
                                           "epsilonAB", "kappaAB",
                                           "epsilon_dispersion", "lambda_dispersion",
                                           "mu");
  }
  if (strcmp(argv[1], "src/DispersionFast.cpp") == 0) {
    DispersionSAFT(R, 0.01, 0.01).create_source(argv[1], "Dispersion", "R",
                                                "epsilon_dispersion", "lambda_dispersion");
  }
  if (strcmp(argv[1], "src/AssociationFast.cpp") == 0) {
    AssociationSAFT(R, 0.01, 0.01, 0.01, 0.01).create_source(argv[1], "Association",
                                                             "R", "epsilonAB", "kappaAB",
                                                             "epsilon_dispersion", "lambda_dispersion");
  }
}
