// -*- mode: C++; -*-

#pragma once

#include "Functionals.h"

Functional SaftFluid(double radius, double epsilon, double kappa,
                     double epsdis, double lambda, double lscale, double mu);
Functional SaftFluid2(double radius, double epsilon, double kappa,
                      double epsdis, double lambda, double lscale, double mu);
Functional EntropySaftFluid2(double radius, double epsilon, double kappa,
                             double epsdis, double lambda, double lscale);

Functional WaterSaft(double R, double epsilon_association, double kappa_association,
                     double epsilon_dispersion, double lambda_dispersion,
                     double length_scaling, double mu);

Functional HardSpheresFast(double radius);
Functional HardSpheresWBFast(double radius);
Functional HardSpheresWBm2(double radius);

Functional HardSpheresNoTensor2(double radius);
Functional ContactAtSphere(double radius);
Functional YuWuCorrelationFast(double radius);

Functional Hydrogen();
