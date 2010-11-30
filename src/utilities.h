// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"
#include "GridDescription.h"

struct LiquidProperties {
  double lengthscale, liquid_density, vapor_density, kT;
  double epsilon_dispersion, lambda_dispersion, epsilonAB, kappaAB;
};

extern LiquidProperties water_prop;
static const double atmospheric_pressure = 101325*3.3989316e-14; // in Hartree/bohr^3

// Ultimately, I'd like surface_tension to be smart about finding a
// converged value.  Or maybe just a second function to do that?
double surface_tension(Minimizer min, Functional f, LiquidProperties prop, bool verbose);

double find_density(Functional f, double kT, double nmin, double nmax);
double pressure(Functional f, double kT, double density);

void equation_of_state(FILE *o, Functional f, double kT, double nmin, double nmax);
void other_equation_of_state(FILE *o, Functional f, double kT, double nmin, double nmax);

long peak_memory();
long current_memory();
void reset_peak_memory();
