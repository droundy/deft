// -*- mode: C++; -*-

#pragma once

#include "Functional.h"

#include "MinimalFunctionals.h"

extern Functional kT;
Functional WithTemperature(const Functional &kTnew, const Functional &f);
extern Functional IdealGasOfVeff;
Functional IdealGas();
Functional EntropyOfIdealGasOfVeff();

Functional HardSphereGas(double radius, double mu);
Functional HardSphereGasRF(double radius, double mu);
Functional HardSpheres(double radius);
Functional HardSpheresRF(double radius);
Functional HardSpheresTarazona(double radius);
Functional HardSpheresWB(double radius);
Functional HardSpheresWBnotensor(double radius);
Functional ChemicalPotential(double chemical_potential);
Functional ExternalPotential(const VectorXd &V);

Functional gHS(Functional n, double radius);
Functional gHScarnahan(Functional n, double radius);
Functional DeltaSAFT(double radius, double epsilon, double kappa,
                     double epsdis, double lambdadis);
Functional Xassociation(double radius, double epsilon, double kappa,
                        double epsdis, double lambdadis);
Functional AssociationSAFT(double radius, double epsilon, double kappa,
                           double epsdis, double lambdadis);
Functional SaftFluidSlow(double radius, double epsilon, double kappa,
                         double epsdis, double lambda, double mu);
Functional SaftFluid(double radius, double epsilon, double kappa,
                     double epsdis, double lambda, double mu);
// SaftExcessEnergy is a functional of density, not effective
// potential.
Functional SaftExcessEnergySlow(double R, double epsilon, double kappa,
                                double epsdis, double lambda,
                                double mu);
Functional DispersionSAFTa1(double radius, double epsdis, double lambda);
Functional DispersionSAFTa2(double radius, double epsdis, double lambda);
Functional DispersionSAFT(double radius, double epsdis, double lambda);
Functional SaftEntropy(double R,
		       double epsilon, double kappa,
		       double epsdis, double lambda);

Functional HardSpheresFast(double radius);
Functional HardSpheresRFFast(double radius);
Functional HardSpheresTarazonaFast(double radius);

Functional GaussianPolynomial(double amplitude, double width, int power);

Functional Identity();

Functional EffectivePotentialToDensity();
Functional Gaussian(double width);

Functional StepConvolve(double radius, Expression r = Expression("R"));

Functional ShellConvolve(double radius, Expression r = Expression("R"));

Functional xShellConvolve(double radius, Expression r = Expression("R"));
Functional yShellConvolve(double radius, Expression r = Expression("R"));
Functional zShellConvolve(double radius, Expression r = Expression("R"));

Functional xxShellConvolve(double radius, Expression r = Expression("R"));
Functional yyShellConvolve(double radius, Expression r = Expression("R"));
Functional zzShellConvolve(double radius, Expression r = Expression("R"));
Functional xyShellConvolve(double radius, Expression r = Expression("R"));
Functional yzShellConvolve(double radius, Expression r = Expression("R"));
Functional zxShellConvolve(double radius, Expression r = Expression("R"));

Functional Pow(int power);
Functional PowAndHalf(int powerMinusHalf);
