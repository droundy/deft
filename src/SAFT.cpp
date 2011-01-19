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
                  (0.5*invdiff*R*n2)*zeta*(Functional(1) +
                                           (0.5*invdiff*R*n2)*(1.0/18)));
}

Functional da1_deta(double radius, double epsdis, double lambdainput);
Functional da1_dlam(double radius, double epsdis, double lambdainput);

Functional gSW(double temp, double R, double epsdis0, double lambda) {
  // This is the approximate *contact* density of a square-well fluid.
  // The formula for this is:
  //      gSW = gHS + 0.25/kT*(da1_deta - lambda/(3*eta)*da1_dlambda)

  // First let's give names to a few constants...
  Functional lam = Functional(lambda).set_name("lambda_dispersion");
  Functional epsdis = Functional(epsdis0).set_name("epsilon_dispersion");
  Functional kT = Functional(temp).set_name("kT");
  Functional n3 = StepConvolve(R);

  Functional da1deta = da1_deta(R, epsdis0, lambda);
  Functional da1dlam = da1_dlam(R, epsdis0, lambda);

  //Functional ghs = gHS(n3, R);
  Functional ghs = gHScarnahan(n3, R);
  return ghs + (Functional(0.25)/kT)*(da1deta - lam/(3*n3)*da1dlam);
}

Functional gHScarnahan_simple(Functional n3) {
  // n3 is the "packing fraction" convolved functional.  It may be an
  // "effective packing fraction", in the case of SAFT-VR.
  return (1 - 0.5*n3)/Pow(3)(1 - n3);
}

Functional gHScarnahan(Functional n3, double R) {
  // n3 is the "packing fraction" convolved functional.  It may be an
  // "effective packing fraction", in the case of SAFT-VR.
  Functional zeta = getzeta(R);
  // The zeta I'm adding here is pretty arbitrary... :(
  return (1 - 0.5*n3)/Pow(3)(1 - n3);
  return zeta*(1 - 0.5*n3)/Pow(3)(1 - n3);
}

Functional dgHScarnahan_dn(Functional n3, double R) {
  // n3 is the "packing fraction" convolved functional.  It may be an
  // "effective packing fraction", in the case of SAFT-VR.
  Functional zeta = getzeta(R);
  // The zeta I'm adding here is pretty arbitrary... :(
  return (2.5-n3)/Pow(4)(1-n3);
  return zeta*(2.5-n3)/Pow(4)(1-n3);
}

Functional DeltaSAFT(double radius, double temperature, double epsilon, double kappa,
                     double epsdis, double lambdadis) {
  Functional g = gSW(temperature, radius, epsdis, lambdadis);
  Functional R(radius);
  R.set_name("R");
  Functional T(temperature);
  T.set_name("kT");
  Functional eps(epsilon);
  eps.set_name("epsilonAB");
  Functional K(kappa);
  K.set_name("kappaAB");
  Functional delta = g*(exp(eps/T) - 1)*K;
  delta.set_name("delta");
  return delta;
}

Functional Xassociation(double radius, double temperature, double epsilon, double kappa,
                        double epsdis, double lambdadis) {
  Functional R(radius);
  R.set_name("R");
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  Functional delta = DeltaSAFT(radius, temperature, epsilon, kappa, epsdis, lambdadis);

  Functional zeta = getzeta(radius);
  Functional X = (sqrt(Functional(1) + 8*n0*zeta*delta) - 1) / (4* n0 * zeta*delta);
  X.set_name("X");
  return X;
}

Functional AssociationSAFT(double radius, double temperature, double epsilon, double kappa,
                           double epsdis, double lambdadis) {
  Functional R(radius);
  R.set_name("R");
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  Functional T(temperature);
  T.set_name("kT");
  Functional zeta = getzeta(radius);
  Functional X = Xassociation(radius, temperature, epsilon, kappa, epsdis, lambdadis);
  return (T*4*n0*zeta*(Functional(0.5) - 0.5*X + log(X))).set_name("association");
}

Functional HardSphereCompressibility(double radius, double temperature) {
  Functional R(radius);
  R.set_name("R");
  Functional kT(temperature);
  kT.set_name("kT");
  return kT;
}

Functional eta_effective(Functional eta, double lambdainput) {
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // The following constants are from equation 37 in Gil-Villegas 1997 paper.
  Functional c1 = Functional(2.25855) - 1.50349*lambda + 0.249434*lambda*lambda;
  Functional c2 = Functional(-0.669270) + 1.40049*lambda - 0.827739*lambda*lambda;
  Functional c3 = Functional(10.1576) - 15.0427*lambda + 5.30827*lambda*lambda;
  // The following equation is equation 36 in Gil-Villegas 1997 paper.
  return c1*eta + c2*sqr(eta) + c3*Pow(3)(eta);
}

Functional detaeff_dlam(Functional eta, double lambdainput) {
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // The following constants are from equation 37 in Gil-Villegas 1997 paper.
  Functional c1 = Functional(-1.50349) + 0.249434*2*lambda;
  Functional c2 = Functional(1.40049) - 0.827739*2*lambda;
  Functional c3 = Functional(-15.0427) + 5.30827*2*lambda;
  // The following equation is equation 36 in Gil-Villegas 1997 paper.
  return c1*eta + c2*sqr(eta) + c3*Pow(3)(eta);
}

Functional detaeff_deta(Functional eta, double lambdainput) {
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // The following constants are from equation 37 in Gil-Villegas 1997 paper.
  Functional c1 = Functional(2.25855) - 1.50349*lambda + 0.249434*lambda*lambda;
  Functional c2 = Functional(-0.669270) + 1.40049*lambda - 0.827739*lambda*lambda;
  Functional c3 = Functional(10.1576) - 15.0427*lambda + 5.30827*lambda*lambda;
  // The following equation is equation 36 in Gil-Villegas 1997 paper.
  return c1 + 2*c2*eta + 3*c3*sqr(eta);
}

Functional eta_for_dispersion(double radius, double lambdainput) {
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  Expression lambdaE("lambda_dispersion");
  lambdaE.set_type("double");
  Expression R("R");
  R.set_type("double");

  return
    (StepConvolve(2*lambdainput*radius, 2*lambdaE*R)
     /(8*Pow(3)(lambda))).set_name("eta_dispersion");

  /*
  return
    ((StepConvolve(2*lambdainput*radius, 2*lambdaE*R) - StepConvolve(radius))
     /(8*Pow(3)(lambda) - 1)).set_name("eta_dispersion");
  */
  
  /*
  return (StepConvolve(lambdainput*radius,lambdaE*R)/Pow(3)(lambda)).set_name("eta_dispersion");
  // FIXME: I think maybe I actually want to compute eta with a larger
  // radius, so as to effectively give the interaction a larger
  // radius? Maybe lambda*radius?
  */
}

Functional DispersionSAFTa1(double radius, double epsdis, double lambdainput) {
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // In Gil-Villegas 1997 paper, packing fraction is called eta...
  Functional eta = eta_for_dispersion(radius, lambdainput);
  eta.set_name("eta");
  Functional eta_eff = eta_effective(eta, lambdainput);
  Functional epsilon_dispersion(epsdis);
  epsilon_dispersion.set_name("epsilon_dispersion");
  // The following equation is equation 35 in Gil-Villegas 1997 paper.
  Functional a1vdw = -4*(Pow(3)(lambda) - 1)*epsilon_dispersion*eta;
  // The following equation is equation 34 in Gil-Villegas 1997 paper.
  return (a1vdw*gHScarnahan(eta_eff, radius)).set_name("a1");
}

Functional da1_dlam(double radius, double epsdis, double lambdainput) {
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // In Gil-Villegas 1997 paper, packing fraction is called eta...
  Functional eta = eta_for_dispersion(radius, lambdainput);
  eta.set_name("eta");
  Functional eta_eff = eta_effective(eta, lambdainput);
  Functional epsilon_dispersion(epsdis);
  epsilon_dispersion.set_name("epsilon_dispersion");
  // The following equation is equation 35 in Gil-Villegas 1997 paper.
  Functional a1vdw_nolam = -4*epsilon_dispersion*eta;
  // The following equation is equation 34 in Gil-Villegas 1997 paper.
  return a1vdw_nolam*(3*sqr(lambda)*gHScarnahan(eta_eff, radius) +
                      (Pow(3)(lambda) - 1)*dgHScarnahan_dn(eta_eff, radius)*
                      detaeff_dlam(eta, lambdainput));
  
}

Functional da1_deta(double radius, double epsdis, double lambdainput) {
  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  // In Gil-Villegas 1997 paper, packing fraction is called eta...
  Functional eta = eta_for_dispersion(radius, lambdainput);
  eta.set_name("eta");
  Functional eta_eff = eta_effective(eta, lambdainput);
  Functional epsilon_dispersion(epsdis);
  epsilon_dispersion.set_name("epsilon_dispersion");
  // The following equation is equation 35 in Gil-Villegas 1997 paper.
  Functional a1vdw_over_eta = -4*(Pow(3)(lambda) - 1)*epsilon_dispersion;
  // The following equation is equation 34 in Gil-Villegas 1997 paper.
  return a1vdw_over_eta*(gHScarnahan(eta_eff, radius) +
                         eta*dgHScarnahan_dn(eta_eff, radius)*
                         detaeff_deta(eta, lambdainput));
  
}

Functional DispersionSAFTa2(double radius, double epsdis, double lambdainput) {
  // In Gil-Villegas 1997 paper, packing fraction is called eta...
  Functional eta = eta_for_dispersion(radius, lambdainput);
  Functional simple_eta_effective = eta_effective(Identity(), lambdainput);
  // The following equation is equation 35 in Gil-Villegas 1997 paper.
  // Actually, it's slightly modified, since the n0 below cancels out
  // the packing fraction by giving us a per-volume rather than
  // per-monomer energy.
  Functional epsilon_dispersion(epsdis);
  epsilon_dispersion.set_name("epsilon_dispersion");
  Functional a1prime = da1_deta(radius, epsdis, lambdainput);

  Functional one_minus_eta = Functional(1) - eta;
  // The following is the Percus-Yevick hard-sphere compressibility
  // factor, see Equation 16 in Gloor 2004 paper.
  //k_hs=(1.d0-eta)**4.0d0/(1.0d0+4.0d0*(eta+eta**2));
  Functional Khs = Pow(4)(one_minus_eta)/(Functional(1) + 4*(eta + sqr(eta)));
  return 0.5*epsilon_dispersion*Khs*eta*a1prime;
}


Functional DispersionSAFT(double radius, double temperature, double epsdis, double lambdainput) {
  Functional kT(temperature);
  kT.set_name("kT");

  Functional lambda(lambdainput);
  lambda.set_name("lambda_dispersion");
  Expression lambdaE("lambda_dispersion");
  lambdaE.set_type("double");
  Expression RE("R");
  RE.set_type("double");
  Functional R(radius);
  R.set_name("R");
  // ndisp is the density of molecules that are within dispersion
  // interaction range of this point.

  /*
  Functional ndisp =
    (StepConvolve(lambdainput*radius,lambdaE*RE) - StepConvolve(radius))
    / ((Pow(3)(lambda) - 1)*((4*M_PI/3)*Pow(3)(R)));
  */
  /*
  Functional ndisp =
    StepConvolve(lambdainput*radius,lambdaE*RE)/(Pow(3)(lambda)*(4*M_PI/3)*Pow(3)(R));
  */
  Functional ndisp = Identity();

  /*
  Functional R(radius);
  R.set_name("R");
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  // I chose to use n0 here because it's the density that measures
  // "how many monomers are neighboring this point, which is what is
  // more relevant in working out interactions.
  */
  Functional a1 = DispersionSAFTa1(radius, epsdis, lambdainput);
  Functional a2 = DispersionSAFTa2(radius, epsdis, lambdainput);
  return (ndisp*(a1 + a2/kT)).set_name("dispersion");
}

Functional SaftFluidSlow(double R, double kT,
                         double epsilon, double kappa,
                         double epsdis, double lambda,
                         double mu
                         ) {
  Functional n = EffectivePotentialToDensity(kT);
  return HardSpheresWBnotensor(R, kT)(n) + IdealGasOfVeff(kT) +
    ChemicalPotential(mu)(n) +
    AssociationSAFT(R, kT, epsilon, kappa, epsdis, lambda)(n) +
    DispersionSAFT(R, kT, epsdis, lambda)(n);
}
