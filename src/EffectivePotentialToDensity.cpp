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
#include "Testables.h"

Functional EffectivePotentialToDensity() {
  Functional Veff = Identity().set_name("Veff");
  return exp(Veff/(-kT()));
}


// The following is a utility function used to create debug prints on
// functionals.
class OfEffectivePotentialClass : public FunctionalInterface {
public:
  OfEffectivePotentialClass(const Functional &myf) : f(myf) {};
  bool append_to_name(const std::string) {
    return false;
  }
  bool I_have_analytic_grad() const {
    return f.I_have_analytic_grad();
  }

  VectorXd transform(const GridDescription &gd, double kT, const VectorXd &data) const {
    return f(gd, kT, (data/(-kT)).cwise().exp());
  }
  double integral(const GridDescription &gd, double kT, const VectorXd &x) const {
    return f.integral(gd, kT, (x/(-kT)).cwise().exp());
  }
  double transform(double kT, double n) const {
    return f(kT, exp(-n/kT));
  }
  double derive(double kT, double V) const {
    double n = exp(-V/kT);
    return (n/-kT)*f.derive(kT, n);
  }
  double d_by_dT(double kT, double n) const {
    return f.d_by_dT(kT, exp(-n/kT));
  }
  Functional grad(const Functional &ingrad, const Functional &V, bool ispgrad) const {
    Functional n = exp(-V/kT());
    if (ispgrad) return f.grad(-ingrad/kT(), n, false);
    return (-n/kT())*f.grad(ingrad, n, false);
  }
  Functional grad_T(const Functional &ingradT) const {
    return OfEffectivePotential(f.grad_T(ingradT));
  }
  void grad(const GridDescription &gd, double kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    if (outpgrad) {
      Grid g(gd);
      g.setZero();
      f.grad(gd, kT, (data/(-kT)).cwise().exp(), ingrad, &g, 0);
      *outgrad += (data/(-kT)).cwise().exp().cwise()*g/(-kT);
      *outpgrad += g/(-kT);
    } else {
      Grid g(gd);
      g.setZero();
      f.grad(gd, kT, (data/(-kT)).cwise().exp(), ingrad, &g, 0);
      *outgrad += (data/(-kT)).cwise().exp().cwise()*g/(-kT);
    }
  }
  void print_summary(const char *prefix, double e, std::string name) const {
    f.print_summary(prefix, e, name);
  }
private:
  Functional f;
  const std::string pattern, args;
};


Functional OfEffectivePotential(const Functional &f) {
  return Functional(new OfEffectivePotentialClass(f));
}
