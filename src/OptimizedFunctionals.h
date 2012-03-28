// -*- mode: C++; -*-

#pragma once

#include "Functionals.h"

Functional SaftFluid(double radius, double epsilon, double kappa,
                     double epsdis, double lambda, double lscale, double mu);
Functional SaftFluid2(double radius, double epsilon, double kappa,
                      double epsdis, double lambda, double lscale, double mu);

Functional Association2(double radius, double epsilon, double kappa,
                        double epsdis, double lambda, double lscale);
Functional Dispersion2(double radius, double epsdis, double lambda, double lscale);

Functional HardSpheresFast(double radius);
Functional HardSpheresRFFast(double radius);
Functional HardSpheresTarazonaFast(double radius);
Functional HardSpheresWBFast(double radius);
Functional HardSpheresWBm2(double radius);

Functional HardSpheresNoTensor2(double radius);
Functional ContactAtSphere(double radius);
Functional YuWuCorrelationFast(double radius);

Functional Hydrogen();