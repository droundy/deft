// -*- mode: C++; -*-

#pragma once

#include "Functional.h"

extern Functional kT;

// HardSpheresNoTensor is used in various SAFT functionals that are
// themselves optimized.
Functional HardSpheresNoTensor(double radius);
  Functional HardSpheresNoTensorGrad(double radius);
  Functional HardSpheresNoTensor_by_dT(double radius);

// Association is the slowest portion of the SAFT free energy, so we
// *really* want to compile it separately.  It'd be beautiful to break
// it into smaller chunks itself, but that looks challenging.
Functional Association(double radius, double epsilon, double kappa,
                       double epsdis, double lambdadis);
  Functional AssociationGrad(double radius, double epsilon, double kappa,
                             double epsdis, double lambdadis);
  Functional Association_by_dT(double radius, double epsilon, double kappa,
                               double epsdis, double lambdadis);

// Dispersion is relatively cheap, but we may as well optimize it
// separately, as it is reused in different places.  It has work in
// common with Association, but I'd pay quite a runtime penalty if it
// could make Association cheaper to compile.
Functional Dispersion(double radius, double epsdis, double lambda);
  Functional DispersionGrad(double radius, double epsdis, double lambda);
  Functional Dispersion_by_dT(double radius, double epsdis, double lambda);

// SaftExcessEnergy is used by the crossover functional, which we'd
// like to be able to optimize.
Functional SaftExcessEnergy(double R, double epsilon, double kappa,
                            double epsdis, double lambda,
                            double mu);
  Functional SaftExcessEnergyGrad(double R, double epsilon, double kappa,
                                  double epsdis, double lambda,
                                  double mu);
  Functional SaftExcessEnergy_by_dT(double R, double epsilon, double kappa,
                                    double epsdis, double lambda,
                                    double mu);
