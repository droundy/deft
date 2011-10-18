// -*- mode: C++; -*-

#pragma once

#include "Functionals.h"

Functional SaftFluid(double radius, double epsilon, double kappa,
                     double epsdis, double lambda, double lscale, double mu);

Functional HardSpheresFast(double radius);
Functional HardSpheresRFFast(double radius);
Functional HardSpheresTarazonaFast(double radius);
Functional HardSpheresWBFast(double radius);
