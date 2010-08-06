// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"
#include "FieldFunctional.h"

Functional IdealGas(const GridDescription &g, double temperature);
Functional ChemicalPotential(const GridDescription &g, double chemical_potential);
Functional ExternalPotential(const Grid &V);

Functional GaussianPolynomial(const GridDescription &g, double amplitude, double width, int power);

FieldFunctional EffectivePotentialToDensity(double temperature);
