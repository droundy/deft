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

Functional HardSphereCompressibility(double radius, double temperature) {
  Functional R(radius);
  R.set_name("R");
  Functional kT(temperature);
  kT.set_name("kT");
  return kT;
}

Functional DispersionSAFTa1(double radius, double temperature, double epsdis, double lambdainput) {
  Functional R(radius);
  R.set_name("R");
  Functional kT(temperature);
  kT.set_name("kT");
  Functional n3 = StepConvolve(radius);
  // FIXME: I think maybe I actually want to compute eta with a larger
  // radius, so as to effectively give the interaction a larger
  // radius? Maybe lambda*radius?
  Functional eta = n3; // In Gil-Villegas 1997 paper, packing fraction is called eta...
  eta.set_name("eta");
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // The following constants are from equation 37 in Gil-Villegas 1997 paper.
  Functional c1 = Functional(2.25855) - 1.50349*lambda + 0.249434*lambda*lambda;
  Functional c2 = Functional(-0.669270) + 1.40049*lambda - 0.827739*lambda*lambda;
  Functional c3 = Functional(10.1576) - 15.0427*lambda + 5.30827*lambda*lambda;
  // The following equation is equation 36 in Gil-Villegas 1997 paper.
  Functional eta_effective = c1*eta + c2*sqr(eta) + c3*Pow(3)(eta);
  eta_effective.set_name("eta_effective");
  // The following equation is equation 35 in Gil-Villegas 1997 paper.
  // Actually, it's slightly modified, since the n0 below cancels out
  // the packing fraction by giving us a per-volume rather than
  // per-monomer energy.
  Functional epsilon_dispersion(epsdis);
  epsilon_dispersion.set_name("epsilon_dispersion");
  Functional a1vdw = -4*(lambda*lambda*lambda - 1)*epsilon_dispersion;
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  // I chose to use n0 here because it's the density that measures
  // "how many monomers are neighboring this point, which is what is
  // more relevant in working out interactions.
  Functional a1 = (n0*a1vdw*gHS(eta_effective, radius)).set_name("a1");
  a1 = (gHScarnahan(eta_effective, radius)).set_name("a1");
  return a1.set_name("A1");
}

Functional DispersionSAFTa2(double radius, double temperature, double epsdis, double lambdainput) {
  Functional R(radius);
  R.set_name("R");
  Functional kT(temperature);
  kT.set_name("kT");
  Functional n3 = StepConvolve(radius);
  // FIXME: I think maybe I actually want to compute eta with a larger
  // radius, so as to effectively give the interaction a larger
  // radius? Maybe lambda*radius?
  Functional eta = n3; // In Gil-Villegas 1997 paper, packing fraction is called eta...
  eta.set_name("eta");
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // The following constants are from equation 37 in Gil-Villegas 1997 paper.
  Functional c1 = Functional(2.25855) - 1.50349*lambda + 0.249434*lambda*lambda;
  Functional c2 = Functional(-0.669270) + 1.40049*lambda - 0.827739*lambda*lambda;
  Functional c3 = Functional(10.1576) - 15.0427*lambda + 5.30827*lambda*lambda;
  // The following equation is equation 36 in Gil-Villegas 1997 paper.
  Functional eta_effective = c1*eta + c2*sqr(eta) + c3*Pow(3)(eta);
  eta_effective.set_name("eta_effective");
  Functional simple_eta_effective = c1*Identity() + c2*Pow(2) + c3*Pow(3);
  // The following equation is equation 35 in Gil-Villegas 1997 paper.
  // Actually, it's slightly modified, since the n0 below cancels out
  // the packing fraction by giving us a per-volume rather than
  // per-monomer energy.
  Functional epsilon_dispersion(epsdis);
  epsilon_dispersion.set_name("epsilon_dispersion");
  Functional a1vdw = -4*(lambda*lambda*lambda - 1)*epsilon_dispersion;
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  // I chose to use n0 here because it's the density that measures
  // "how many monomers are neighboring this point, which is what is
  // more relevant in working out interactions.
  Functional a1 = (n0*a1vdw*gHScarnahan(eta_effective, radius)).set_name("a1");
  Functional gHSprime = gHScarnahan(simple_eta_effective, radius).grad(Identity(), eta, false);
  // FIXME: for some reason, setting the name of gHSprime causes a problem...  :(
  //gHSprime.set_name("gHSprime");
  Functional one_minus_eta = Functional(1) - eta;
  // The following is the Percus-Yevick hard-sphere compressibility
  // factor, see Equation 16 in Gloor 2004 paper.
  Functional Khs = Pow(4)(one_minus_eta)/(Functional(1) + 4*eta + 4*sqr(eta));
  // a2 is 1/2 beta epsilon Khs \frac{\partial a_1}{\partial \eta}
  Functional a2 = 0.5*epsilon_dispersion*Khs*
    (a1 - 4*(lambda*lambda*lambda-1)*epsilon_dispersion*n0*sqr(eta)*gHSprime);
  return (a2/kT).set_name("A2");
}

Functional DispersionSAFT(double radius, double temperature, double epsdis, double lambdainput) {
  Functional R(radius);
  R.set_name("R");
  Functional kT(temperature);
  kT.set_name("kT");
  Functional n3 = StepConvolve(radius);
  // FIXME: I think maybe I actually want to compute eta with a larger
  // radius, so as to effectively give the interaction a larger
  // radius? Maybe lambda*radius?
  Functional eta = n3; // In Gil-Villegas 1997 paper, packing fraction is called eta...
  eta.set_name("eta");
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // The following constants are from equation 37 in Gil-Villegas 1997 paper.
  Functional c1 = Functional(2.25855) - 1.50349*lambda + 0.249434*lambda*lambda;
  Functional c2 = Functional(-0.669270) + 1.40049*lambda - 0.827739*lambda*lambda;
  Functional c3 = Functional(10.1576) - 15.0427*lambda + 5.30827*lambda*lambda;
  // The following equation is equation 36 in Gil-Villegas 1997 paper.
  Functional eta_effective = c1*eta + c2*sqr(eta) + c3*Pow(3)(eta);
  eta_effective.set_name("eta_effective");
  Functional simple_eta_effective = c1*Identity() + c2*Pow(2) + c3*Pow(3);
  // The following equation is equation 35 in Gil-Villegas 1997 paper.
  // Actually, it's slightly modified, since the n0 below cancels out
  // the packing fraction by giving us a per-volume rather than
  // per-monomer energy.
  Functional epsilon_dispersion(epsdis);
  epsilon_dispersion.set_name("epsilon_dispersion");
  Functional a1vdw = -4*(lambda*lambda*lambda - 1)*epsilon_dispersion;
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  // I chose to use n0 here because it's the density that measures
  // "how many monomers are neighboring this point, which is what is
  // more relevant in working out interactions.
  Functional a1 = (a1vdw*gHScarnahan(eta_effective, radius)).set_name("a1");
  Functional gHSprime = gHScarnahan(simple_eta_effective, radius).grad(Identity(), eta, false);
  // FIXME: for some reason, setting the name of gHSprime causes a problem...  :(
  //gHSprime.set_name("gHSprime");
  Functional one_minus_eta = Functional(1) - eta;
  // The following is the Percus-Yevick hard-sphere compressibility
  // factor, see Equation 16 in Gloor 2004 paper.
  Functional Khs = Pow(4)(one_minus_eta)/(Functional(1) + 4*eta + 4*sqr(eta));
  // a2 is 1/2 beta epsilon Khs \frac{\partial a_1}{\partial \eta}
  Functional a2 = 0.5*epsilon_dispersion*Khs*
    (a1 + 4*(1 - lambda*lambda*lambda)*epsilon_dispersion*sqr(eta)*gHSprime);
  return (n0*(a1 + a2/kT)).set_name("dispersion");
}

Functional SaftFluidSlow(double R, double kT,
                         double epsilon, double kappa,
                         double epsdis, double lambda,
                         double mu
                         ) {
  Functional n = EffectivePotentialToDensity(kT);
  return HardSpheresWBnotensor(R, kT)(n) + IdealGasOfVeff(kT) +
    ChemicalPotential(mu)(n) +
    AssociationSAFT(R, kT, epsilon, kappa)(n) +
    DispersionSAFT(R, kT, epsdis, lambda)(n);
}
