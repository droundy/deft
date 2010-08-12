// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"
#include "FieldFunctional.h"

Functional IdealGas(double temperature);
Functional ChemicalPotential(double chemical_potential);
Functional ExternalPotential(const Grid &V);

Functional GaussianPolynomial(double amplitude, double width, int power);

FieldFunctional EffectivePotentialToDensity(double temperature);
