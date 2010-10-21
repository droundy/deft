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

const double kT = water_prop.kT; // room temperature in Hartree
const double R = 2.7;
const double mu = 0;

int main(int, char **argv) {
  if (strcmp(argv[1], "src/IdealGasFast.cpp") == 0)
    IdealGas(kT).create_source("src/IdealGasFast.cpp", "IdealGasFast", "kT");
  if (strcmp(argv[1], "src/HardSpheresFast.cpp") == 0)
    HardSpheres(kT, R).create_source("src/HardSpheresFast.cpp", "HardSpheresFast", "R", "kT");
  if (strcmp(argv[1], "src/HardSpheresRFFast.cpp") == 0)
    HardSpheresRF(kT, R).create_source("src/HardSpheresRFFast.cpp", "HardSpheresRFFast", "R", "kT");
  if (strcmp(argv[1], "src/HardSpheresTarazonaFast.cpp") == 0)
    HardSpheresTarazona(kT, R).create_source("src/HardSpheresTarazonaFast.cpp", "HardSpheresTarazonaFast", "R", "kT");
  if (strcmp(argv[1], "src/HardSpheresNoTensor.cpp") == 0)
    HardSpheresWBnotensor(kT, R).create_source("src/HardSpheresNoTensor.cpp", "HardSpheresNoTensor", "R", "kT");
  if (strcmp(argv[1], "src/HardSphereGasFast.cpp") == 0) {
    Functional n = EffectivePotentialToDensity(kT);
    Functional f = HardSpheresRF(R, kT)(n) + IdealGasOfVeff(kT) + ChemicalPotential(mu)(n);
    f.create_source("src/HardSphereGasRFFast.cpp", "HardSphereGasRF", "R", "kT", "mu");

    f = HardSpheres(R, kT)(n) + IdealGasOfVeff(kT) + ChemicalPotential(mu)(n);
    f.create_source("src/HardSphereGasFast.cpp", "HardSphereGas", "R", "kT", "mu");
  }
}
