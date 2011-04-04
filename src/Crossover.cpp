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
#include "Crossover.h"

extern Functional kT;

Functional Crossover(Functional a_res, double Gi_in,
                     double Tc_in, double nc_in,
                     double T0c_in, double n0c_in) {
  //a_res = ReportOn(a_res, "a_res");
  // The notation inside this function comes from kiselev2006new.
  Functional n = Identity();
  Functional nc(nc_in, "nc");
  Functional n0c(n0c_in, "n0c");
  Functional Tc(Tc_in, "Tc");
  //Tc = ReportOn(Tc, "Tc");
  Functional T0c(T0c_in, "T0c");

  // The following definitions follow kiselev2006new equation 5:
  Functional tau = (kT/Tc - Functional(1)).set_name("tau"); // the reduced temperature
  //tau = ReportOn(tau, "tau");
  Functional eta = (nc/n - Functional(1)).set_name("eta"); // the reduced molar volume
  //eta = ReportOn(eta, "eta");
  Functional DeltaTc = (Tc/T0c - Functional(1)).set_name("DeltaTc"); // dimensionless shift of critical temperature
  Functional Delta_vc = (n0c/nc - Functional(1)).set_name("Delta_vc"); // dimensionless shift of molar volume

  Functional v1(1.0, "v1");
  Functional d1(1.0, "d1");
  Functional etahat = eta*(Functional(1) + v1*exp(-10*eta)+d1*tau); // kiselev2006new right after eqn 12
  //etahat = ReportOn(etahat, "etahat");

  const double alpha = 0.11; // A universal exponent
  const double beta = 0.325; // A universal exponent
  const double gamma = 2 - 2*beta - alpha; // A universal exponent
  const double log_one_over_beta = log(1/beta);
  const double Delta1 = 0.51; // A universal exponent
  const double log_2_Delta1 = log(2*Delta1);

  const double b = sqrt(1.359); // The LM parameter from kiselev2006new
  const double m0 = 2.5; // from kiselev2006new before equation 19.
  const double b_over_m0 = b/m0;
  const double b_over_m0_to_oo_beta = pow(b_over_m0, 1/beta);

  Functional etahat_to_oo_beta = exp(log_one_over_beta*log(abs(etahat)));
  //etahat_to_oo_beta = ReportOn(etahat_to_oo_beta, "etahat_to_oo_beta");
  Functional four_b_etc_plus_2_tau = 4*b_over_m0_to_oo_beta*etahat_to_oo_beta + 2*tau;
  four_b_etc_plus_2_tau.set_name("four_b_etc_plus_2_tau");
  //four_b_etc_plus_2_tau = ReportOn(four_b_etc_plus_2_tau, "four_b_etc_plus_2_tau");
  Functional Gi(Gi_in, "Gi");
  // kiselev2006new equation 7:
  Functional q2 =
    (four_b_etc_plus_2_tau + sqrt(sqr(four_b_etc_plus_2_tau) + 12*sqr(tau))) * (Functional(1)/(6*Gi));
  q2.set_name("q_sqr");
  Functional q = sqrt(q2);
  //Functional q = sqrt(q2);

  // kiselev2006new equation 6:
  Functional Y = exp(log_2_Delta1*q/(q+Functional(1.0)));

  // kiselev2006new equation 4:
  Functional tau_bar = tau*exp(-Y*Functional(log(alpha/Delta1)))
    + (tau + Functional(1.0))*DeltaTc*exp(Y*Functional(log(2*(2-alpha)/(3*Delta1))));
  // kiselev2006new equation 5:
  Functional eta_bar = eta*exp(Y*Functional(log((gamma-2*beta)/(4*Delta1))))
    + (eta + Functional(1.0))*Delta_vc*exp(Y*Functional(log((2-alpha)/(2*Delta1))));

  // Here we compute an effective density and effective temperature
  // using the crossover approach.
  Functional n_effective = n0c*(eta_bar + Functional(1.0));
  Functional T_effective = T0c*(tau_bar + Functional(1.0));
  
  // p0 is the pressure of the "classical" system
  Functional p0 = n*kT - n*a_res.grad(Identity(), n, false);
  // p0 is the "dimensionless" pressure scaled according to the critical point.
  Functional p0bar = p0*n0c/kT;

  // In kiselev2006new, \Delta v really ought to be called \eta_0...  :(
  Functional Delta_v = (n0c/n - Functional(1.0)).set_name("Delta_v");

  // kiselev2006new equation 3
  Functional a_res_of_n0c = a_res;
  a_res_of_n0c = a_res_of_n0c.append_to_name("_of_n0c")(n0c); // This is awkward.
  Functional a_bg = a_res_of_n0c + p0*Delta_v; // This omits the ideal gas free energy!
  printf("I am about to create Delta_a\n");
  // kiselev2006new equation 2
  Functional Delta_a = a_res
    - a_res_of_n0c
    + n*p0*Delta_v - n*log(Delta_v + Functional(1));
  // Here we substitute eta_bar and tau_bar into equation 2, as instructed.
  //printf("I am about to RenameAppending and WithTemperature\n");
  //Delta_a = RenameAppending(WithTemperature(T_effective, Delta_a(n_effective)),
  //                          "_of_Teff").set_name("Delta_a");

  printf("Before I go any further, I'll check that we can compute the temperature derivative...\n");
  Delta_a(n_effective).grad_T(Functional(1));

  printf("I am actually about to WithTemperature...\n");
  Delta_a = WithTemperature(T_effective, Delta_a.append_to_name("_of_neff")(n_effective));
  printf("And now I am actually about to RenameAppending...\n");
  Delta_a = Delta_a.append_to_name("_of_Teff").set_name("Delta_a");

  //Delta_a = ReportOn(Delta_a, "Delta_a");

  return Delta_a + a_bg;
}

// The following is a utility function used to create debug prints on
// functionals.
class ReportOnClass : public FunctionalInterface {
public:
  ReportOnClass(const Functional &myf, const char *myname)
    : f(myf), name(myname) {};
  bool append_to_name(const std::string x) {
    f.append_to_name(x);
    return true;
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &kT, const VectorXd &data) const {
    return f(gd, kT, data);
  }
  double transform(double kT, double n) const {
    double fn = f(kT, n);
    printf("\t%s(%g, %g) = %g\n", name, n, kT, fn);
    return fn;
  }
  double derive(double kT, double n) const {
    return f.derive(kT, n);
  }
  Expression derive_homogeneous(const Expression &kT, const Expression &x) const {
    return f.derive_homogeneous(kT, x);
  }
  double d_by_dT(double kT, double n) const {
    return f.d_by_dT(kT, n);
  }
  Functional grad(const Functional &ingrad, const Functional &n, bool ispgrad) const {
    return f.grad(ingrad, n, ispgrad);
  }
  Functional grad_T(const Functional &ingradT) const {
    return f.grad_T(ingradT);
  }
  void grad(const GridDescription &gd, const VectorXd &kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    f.grad(gd, kT, data, ingrad, outgrad, outpgrad);
  }
  Expression printme(const Expression &kT, const Expression &x) const {
    return f.printme(kT, x);
  }
private:
  Functional f;
  const char *name;
};


Functional ReportOn(Functional f, const char *name) {
  printf("I am going to report on %s...\n", name);
  return Functional(new ReportOnClass(f, name));
}
