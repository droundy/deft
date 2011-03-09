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

  VectorXd transform(double kT, const GridDescription &, const VectorXd &data) const {
    return kT*VectorXd::Ones(data.rows());
  }
  double transform(double kT, double) const {
    return kT;
  }
  double derive(double, double) const {
    return 0;
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return 0;
  }
  void grad(double, const GridDescription &, const VectorXd &, const VectorXd &, VectorXd *, VectorXd *) const {
  }
  Expression printme(const Expression &) const {
    return Expression("kT").set_type("double");
  }
};

Functional kT = Functional(new Temperature(), "kT");


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

static Functional CreateIdealGasOfVeff() {
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
