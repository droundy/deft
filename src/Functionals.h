// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"

Functional IdealGas(double temperature);
Functional IdealGasOfVeff(double temperature);
Functional HardSphereGas(double radius, double temperature, double mu);
Functional HardSpheres(double radius, double temperature);
Functional HardSpheresRF(double radius, double temperature);
Functional HardSpheresTarazona(double radius, double temperature);
Functional HardSpheresWB(double radius, double temperature);
Functional HardSpheresWBnotensor(double radius, double temperature);
Functional ChemicalPotential(double chemical_potential);
Functional ExternalPotential(const VectorXd &V);

Functional HardSpheresFast(double radius, double temperature);
Functional HardSpheresRFFast(double radius, double temperature);
Functional HardSpheresTarazonaFast(double radius, double temperature);
Functional HardSpheresNoTensor(double radius, double temperature);

Functional GaussianPolynomial(double amplitude, double width, int power);

Functional Identity();

Functional EffectivePotentialToDensity(double temperature);
Functional Gaussian(double width);

Functional StepConvolve(double radius);

Functional ShellConvolve(double radius);

Functional xShellConvolve(double radius);
Functional yShellConvolve(double radius);
Functional zShellConvolve(double radius);

Functional xxShellConvolve(double radius);
Functional yyShellConvolve(double radius);
Functional zzShellConvolve(double radius);
Functional xyShellConvolve(double radius);
Functional yzShellConvolve(double radius);
Functional zxShellConvolve(double radius);

Functional Pow(int power);
