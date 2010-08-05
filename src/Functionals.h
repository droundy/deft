// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"
#include "FieldFunctional.h"

Functional IdealGas(const GridDescription &g, double temperature);
Functional ChemicalPotential(const GridDescription &g, double chemical_potential);
Functional ExternalPotential(const Grid &V);

FieldFunctional EffectivePotentialToDensity(double temperature);
