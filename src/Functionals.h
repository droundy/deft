// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"

Functional IdealGas(double temperature);
Functional IdealGasOfVeff(double temperature);
Functional EntropyOfIdealGasOfVeff(double T);
Functional HardSphereGas(double radius, double temperature, double mu);
Functional HardSphereGasRF(double radius, double temperature, double mu);
Functional HardSpheres(double radius, double temperature);
Functional HardSpheresRF(double radius, double temperature);
Functional HardSpheresTarazona(double radius, double temperature);
Functional HardSpheresWB(double radius, double temperature);
Functional HardSpheresWBnotensor(double radius, double temperature);
Functional ChemicalPotential(double chemical_potential);
Functional ExternalPotential(const VectorXd &V);

Functional gHS(Functional n, double radius);
Functional gHScarnahan(Functional n, double radius);
Functional DeltaSAFT(double radius, double temperature, double epsilon, double kappa,
                     double epsdis, double lambdadis);
Functional Xassociation(double radius, double temperature, double epsilon, double kappa,
                        double epsdis, double lambdadis);
Functional AssociationSAFT(double radius, double temperature, double epsilon, double kappa,
                           double epsdis, double lambdadis);
Functional SaftFluidSlow(double radius, double temperature, double epsilon, double kappa,
                         double epsdis, double lambda, double mu);
Functional SaftFluid(double radius, double temperature,
                     double epsilon, double kappa,
                     double epsdis, double lambda, double mu);
Functional DispersionSAFTa1(double radius, double epsdis, double lambda);
Functional DispersionSAFTa2(double radius, double epsdis, double lambda);
Functional DispersionSAFT(double radius, double temperature, double epsdis, double lambda);
Functional SaftEntropy(double R, double temp,
		       double epsilon, double kappa,
		       double epsdis, double lambda);

Functional HardSpheresFast(double radius, double temperature);
Functional HardSpheresRFFast(double radius, double temperature);
Functional HardSpheresTarazonaFast(double radius, double temperature);
Functional HardSpheresNoTensor(double radius, double temperature);

Functional GaussianPolynomial(double amplitude, double width, int power);

Functional Identity();

Functional EffectivePotentialToDensity(double temperature);
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
