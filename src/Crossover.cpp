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

extern Functional kT;

Functional Crossover(Functional fexc0, double Gi_in,
                     double Tc_in, double nc_in,
                     double T0c_in, double n0c_in) {
  // The notation inside this function comes from kiselev2006new.
  Functional n = Identity();
  Functional nc(nc_in, "nc");
  Functional n0c(n0c_in, "n0c");
  Functional Tc(Tc_in, "Tc");
  Functional T0c(T0c_in, "T0c");

  // The following definitions follow kiselev2006new equation 5:
  Functional tau = (kT/Tc - 1).set_name("tau"); // the reduced temperature
  Functional eta = (nc/n - 1).set_name("eta"); // the reduced molar volume
  Functional DeltaTc = (Tc/T0c - 1).set_name("DeltaTc"); // dimensionless shift of critical temperature
  Functional Delta_vc = (n0c/nc - 1).set_name("Delta_vc"); // dimensionless shift of molar volume

  Functional etahat = eta; // FIXME

  const double alpha = 0.11; // A universal exponent
  const double beta = 0.325; // A universal exponent
  const double gamma = 2 - 2*beta - alpha; // A universal exponent
  const double log_one_over_beta = log(1/beta);
  const double Delta1 = 0.51; // A universal exponent
  const double log_2_Delta1 = log(2*Delta1);

  const double b = sqrt(1.359); // The LM parameter from kiselev2006new
  const double m0 = 2.5; // from kiselev2006new before equation 19.
  const double b_over_m0 = b/m0; // FIXME
  const double b_over_m0_to_oo_beta = pow(b_over_m0, 1/beta);

  Functional etahat_to_oo_beta = exp(log_one_over_beta*log(etahat));
  Functional four_b_etc_plus_2_tau = 4*b_over_m0_to_oo_beta*etahat_to_oo_beta + 2*tau;
  four_b_etc_plus_2_tau.set_name("four_b_etc_plus_2_tau");
  Functional Gi(Gi_in, "Gi");
  // kiselev2006new equation 7:
  Functional q2 =
    (four_b_etc_plus_2_tau + sqrt(sqr(four_b_etc_plus_2_tau) + 12*sqr(tau))) * (Functional(1)/(6*Gi));
  q2.set_name("q_sqr");
  Functional q = sqrt(q);

  // kiselev2006new equation 6:
  Functional Y = exp(log_2_Delta1*q/(q+1));

  // kiselev2006new equation 4:
  Functional tau_bar = tau*exp(-Y*log(alpha/Delta1))
    + (tau + 1)*DeltaTc*exp(Y*log(2*(2-alpha)/(3*Delta1)));
  // kiselev2006new equation 5:
  Functional eta_bar = eta*exp(Y*log((gamma-2*beta)/(4*Delta1)))
    + (eta + 1)*Delta_vc*exp(Y*log((2-alpha)/(2*Delta1)));

  // Here we compute an effective density and effective temperature
  // using the crossover approach.
  Functional n_effective = n0c*(eta_bar + 1);
  Functional T_effective = T0c*(tau_bar + 1);
  
  Functional a_res;
  // p0 is the pressure of the "classical" system
  Functional p0 = n*kT - n*fexc0.grad(Identity(), n, false);
  // p0 is the "dimensionless" pressure scaled according to the critical point.
  Functional p0bar = p0*n0c/kT;

  // In kiselev2006new, \Delta v really ought to be called \eta_0...  :(
  Functional Delta_v = (n0c/n - 1).set_name("Delta_v");

  // kiselev2006new equation 3
  Functional a_bg = a_res(n0c) + p0*Delta_v; // This omits the ideal gas free energy!
  // kiselev2006new equation 2
  Functional Delta_a = fexc0 - fexc0(n0c) + p0*Delta_v - log(Delta_v + 1);
  // Here we substitute eta_bar and tau_bar into equation 2, as instructed.
  Delta_a = WithTemperature(T_effective, Delta_a(n_effective)).set_name("Delta_a");

  return Delta_a + a_bg;
}
