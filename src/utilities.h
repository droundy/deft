// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"
#include "GridDescription.h"

struct LiquidProperties {
  double lengthscale, liquid_density, vapor_density, kT;
};

extern LiquidProperties water_prop;

// Ultimately, I'd like surface_tension to be smart about finding a
// converged value.  Or maybe just a second function to do that?
double surface_tension(Minimizer min, Functional f, LiquidProperties prop, bool verbose);

long peak_memory();
