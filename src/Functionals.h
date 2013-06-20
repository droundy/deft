// -*- mode: C++; -*-

#pragma once

#include "Functional.h"

#include "MinimalFunctionals.h"

Functional IdealGasOfVeff();
Functional IdealGas();
Functional EntropyOfIdealGasOfVeff();
Functional EntropyOfIdealGas();

Functional HardSphereGas(double radius, double mu);
Functional HardSphereGasRF(double radius, double mu);
Functional HardSpheres(double radius);
Functional HardSpheresRF(double radius);
Functional HardSpheresTarazona(double radius);
Functional HardSpheresWB(double radius);
Functional HardSpheresWBnotensor(double radius);
Functional HardSpheresWBm2slow(double radius);
Functional ChemicalPotential(double chemical_potential);
Functional ExternalPotential(const VectorXd &V);

Functional gHS(Functional n, double radius);
Functional gHScarnahan(Functional n, double radius);
Functional DeltaSAFT(double radius, double epsilon, double kappa,
                     double epsdis, double lambdadis, double lscale);
Functional Xassociation(double radius, double epsilon, double kappa,
                        double epsdis, double lambdadis, double lscale);
Functional AssociationSAFT(double radius, double epsilon, double kappa,
                           double epsdis, double lambdadis, double lscale);
Functional SaftFluidSlow(double radius, double epsilon, double kappa,
                         double epsdis, double lambda, double lscale, double mu);

Functional DispersionSAFTa1(double radius, double epsdis, double lambda, double lscale);
Functional DispersionSAFTa2(double radius, double epsdis, double lambda, double lscale);
Functional DispersionSAFT(double radius, double epsdis, double lambda, double lscale);
Functional SaftEntropy(double R,
		       double epsilon, double kappa,
		       double epsdis, double lambda, double lscale);

Functional GaussianPolynomial(double amplitude, double width, int power);

Functional Identity();

Functional EffectivePotentialToDensity();
Functional OfEffectivePotential(const Functional &f);
Functional Gaussian(double width);

Functional GaussianConvolve(double width);
Functional StepConvolve(double radius);

Functional ShellConvolve(double radius);
Functional ShellPrimeConvolve(double radius);

Functional xShellConvolve(double radius);
Functional yShellConvolve(double radius);
Functional zShellConvolve(double radius);

Functional xShellPrimeConvolve(double radius);
Functional yShellPrimeConvolve(double radius);
Functional zShellPrimeConvolve(double radius);

Functional xxShellConvolve(double radius);
Functional yyShellConvolve(double radius);
Functional zzShellConvolve(double radius);
Functional xyShellConvolve(double radius);
Functional yzShellConvolve(double radius);
Functional zxShellConvolve(double radius);

Functional Pow(int power);
Functional PowAndHalf(int powerMinusHalf);
