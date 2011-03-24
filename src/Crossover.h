// -*- mode: C++; -*-

#pragma once

#include "Functionals.h"

Functional Crossover(Functional fexcess, double Gi, double Tc, double nc, double T0c, double n0c);

Functional ReportOn(Functional f, const char *name);
