// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Deft Authors
//
// Deft is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with deft; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include "Functionals.h"
#include <stdio.h>
#include <math.h>

static Functional getzeta(double radius) {
  Functional n2 = ShellConvolve(radius);
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  return 1 - (sqr(n2x) + sqr(n2y) + sqr(n2z))/sqr(n2);
  //return (1 - (sqr(n2x) + sqr(n2y) + sqr(n2z))/sqr(n2)).set_name("xi");
}

Functional xi2(double radius) {
  Functional n2 = ShellConvolve(radius);
  Functional R(radius);
  R.set_name("R");
  return (M_PI*sqr(R)*n2).set_name("xi2");
}

Functional xi3(double radius) {
  Functional n3 = StepConvolve(radius);
  Functional R(radius);
  R.set_name("R");
  return ((M_PI*4/3)*Pow(3)(R)*n3).set_name("xi3");
}

Functional gHS(Functional n3, double R) {
  // n3 is the "packing fraction" convolved functional.  It may be an
  // "effective packing fraction", in the case of SAFT-VR.
  Functional zeta = getzeta(R);
  Functional n2 = ShellConvolve(R);
  Functional invdiff = Functional(1)/(1-n3);
  return invdiff*(Functional(1) +
                  invdiff*n2*zeta*(Functional(1.0/4) +
                                   invdiff*(1.0/18/4)*n2));
}

Functional gHScarnahan(Functional n3, double R) {
  // n3 is the "packing fraction" convolved functional.  It may be an
  // "effective packing fraction", in the case of SAFT-VR.
  Functional zeta = getzeta(R);
  // The zeta I'm adding here is pretty arbitrary... :(
  return (1 - 0.5*n3)/Pow(3)(1 - n3);
  return zeta*(1 - 0.5*n3)/Pow(3)(1 - n3);
}

Functional DeltaSAFT(double radius, double temperature, double epsilon, double kappa) {
  Functional n3 = StepConvolve(radius);
  Functional g = gHS(n3, radius);
  Functional R(radius);
  R.set_name("R");
  Functional T(temperature);
  T.set_name("kT");
  Functional eps(epsilon);
  eps.set_name("epsilonAB");
  Functional K(kappa);
  K.set_name("kappaAB");
  Functional delta = g*(exp(eps/T) - 1)*8*Pow(3)(R)*K;
  delta.set_name("delta");
  return delta;
}

Functional Xassociation(double radius, double temperature, double epsilon, double kappa) {
  Functional R(radius);
  R.set_name("R");
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  Functional delta = DeltaSAFT(radius, temperature, epsilon, kappa);

  Functional zeta = getzeta(radius);
  Functional X = (sqrt(Functional(1) + 8*n0*zeta*delta) - 1) / (4* n0 * zeta*delta);
  X.set_name("X");
  return X;
}

Functional AssociationSAFT(double radius, double temperature, double epsilon, double kappa) {
  Functional R(radius);
  R.set_name("R");
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  Functional T(temperature);
  T.set_name("kT");
  Functional zeta = getzeta(radius);
  Functional X = Xassociation(radius, temperature, epsilon, kappa);
  return (T*4*n0*zeta*(Functional(0.5) - 0.5*X + log(X))).set_name("association");
}

Functional SaftFluidSlow(double R, double kT,
                         double epsilon, double kappa, double mu) {
  Functional n = EffectivePotentialToDensity(kT);
  return HardSpheresWBnotensor(R, kT)(n) + IdealGasOfVeff(kT) +
    ChemicalPotential(mu)(n) +
    AssociationSAFT(R, kT, epsilon, kappa)(n);
}
