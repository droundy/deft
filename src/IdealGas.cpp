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

class Temperature : public FunctionalInterface {
public:
  Temperature() {
  }

  VectorXd transform(const GridDescription &, const VectorXd &kT, const VectorXd &) const {
    return kT;
  }
  double transform(double kT, double) const {
    return kT;
  }
  double derive(double, double) const {
    return 0;
  }
  double d_by_dT(double, double) const {
    return 1;
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return 0;
  }
  Functional grad_T(const Functional &ingrad) const {
    return ingrad;
  }
  void grad(const GridDescription &, const VectorXd &, const VectorXd &, const VectorXd &,
            VectorXd *, VectorXd *) const {
  }
  Expression printme(const Expression &) const {
    return Expression("kT");
  }
};

Functional kT = Functional(new Temperature(), "kT");

class WithTemperatureClass : public FunctionalInterface {
public:
  WithTemperatureClass(const Functional &newkT, const Functional &myf) : mykT(newkT), f(myf) {};

  VectorXd transform(const GridDescription &gd, const VectorXd &kT, const VectorXd &data) const {
    return f(gd, mykT(gd, kT, data), data);
  }
  double transform(double kT, double n) const {
    return f(mykT(kT, n), n);
  }
  double derive(double kT, double n) const {
    double kTnew = mykT(kT, n);
    // d/dn = d/dn + d/dT * dTnew/dn
    return f.derive(kTnew, n) + f.d_by_dT(kTnew, n)*mykT.derive(kT, n);
  }
  double d_by_dT(double kT, double n) const {
    double kTnew = mykT(kT, n);
    // d/dn = d/dn + d/dT * dTnew/dn
    return f.d_by_dT(kTnew, n)*mykT.d_by_dT(kT, n);
  }
  Functional grad(const Functional &ingrad, const Functional &n, bool ispgrad) const {
    return WithTemperature(mykT, f.grad(ingrad, n, ispgrad))
      + mykT.grad(f.grad_T(ingrad), n, ispgrad);
  }
  Functional grad_T(const Functional &ingradT) const {
    return mykT.grad_T(WithTemperature(mykT, f.grad_T(ingradT)));
  }
  void grad(const GridDescription &gd, const VectorXd &kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    VectorXd kTnew = mykT(gd, kT, data);
    f.grad(gd, kTnew, data, ingrad, outgrad, outpgrad);
    // Now include gradient of any density-dependence in the effective
    // temperature:
    mykT.grad(gd, kT, data, f.grad_T(1)(gd, kTnew, data).cwise()*ingrad, outgrad, outpgrad);
  }
  Expression printme(const Expression &x) const {
    Expression kTnew = mykT.printme(x);
    Expression fnew = f.printme(x);
    fnew.ReplaceThisSubexpression(Expression("kT").set_type("double"), kTnew);
    return fnew;
  }
private:
  Functional mykT, f;
  };

Functional WithTemperature(const Functional &newkT, const Functional &f) {
  return Functional(new WithTemperatureClass(newkT, f));
}

Functional find_nQ() {
  Functional mass = 18*1822.8885; // FIXME: molecular weight of water
  // Note:  hbar is one in atomic units, yay!
  // nQ = (m*kT/(2 pi hbar^2))^3/2
  return PowAndHalf(3)(mass*kT/(2*M_PI));
}

Functional find_dnQ_dT() {
  Functional mass = 18*1822.8885; // FIXME: molecular weight of water
  // Note:  hbar is one in atomic units, yay!
  // nQ = (m*kT/(2 pi hbar^2))^3/2
  return 1.5*PowAndHalf(3)(mass/(2*M_PI))*sqrt(kT);
}

Functional IdealGas() {
  Functional n = Identity();
  Functional nQ = find_nQ();
  return (kT*n*(log(n/nQ) - 1)).set_name("ideal_gas");
}

Functional CreateIdealGasOfVeff() {
  Functional Veff = Identity().set_name("Veff");
  Functional n = exp(-Veff /kT);
  Functional nQ = find_nQ();
  return (-(Veff + kT*log(nQ) + kT)*n).set_name("ideal_gas");
}

Functional IdealGasOfVeff = CreateIdealGasOfVeff();

Functional EntropyOfIdealGasOfVeff() {
  Functional Veff = Identity().set_name("Veff");
  Functional n = exp(-Veff/kT);
  Functional nQ = find_nQ();
  Functional dnQ_dT = find_dnQ_dT();
  // The following is also known as the Sackur-Tetrode equation
  return ((Veff/kT + log(nQ) + 2.5)*n).set_name("ideal_gas_entropy");
}
