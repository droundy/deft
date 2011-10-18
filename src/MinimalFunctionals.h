// -*- mode: C++; -*-

#pragma once

#include "Functional.h"

// This file is intended for functionals that are themselves
// optimized, but are also needed by other optimized functionals.

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
                       double epsdis, double lambdadis, double lscale);
  Functional AssociationGrad(double radius, double epsilon, double kappa,
                             double epsdis, double lambdadis, double lscale);
  Functional Association_by_dT(double radius, double epsilon, double kappa,
                               double epsdis, double lambdadis, double lscale);

// Dispersion is relatively cheap, but we may as well optimize it
// separately, as it is reused in different places.  It has work in
// common with Association, but I'd pay quite a runtime penalty if it
// could make Association cheaper to compile.
Functional Dispersion(double radius, double epsdis, double lambda, double lscale);
  Functional DispersionGrad(double radius, double epsdis, double lambda, double lscale);
  Functional Dispersion_by_dT(double radius, double epsdis, double lambda, double lscale);
