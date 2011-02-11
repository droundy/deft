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

class Choose : public FunctionalInterface {
public:
  Choose(double c, const Functional &fa, const Functional &fb) : cut(c), flow(fa), fhigh(fb) {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    VectorXd out(data);
    const int N = data.cols()*data.rows();
    for (int i=0; i<N; i++) {
      double x = data[i];
      out[i] = (x<cut) ? flow(x) : fhigh(x);
    }
    return out;
  }
  double transform(double x) const {
    return (x<cut) ? flow(x) : fhigh(x);
  }
  double derive(double x) const {
    return (x<cut) ? flow.derive(x) : fhigh.derive(x);
  }

  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return ingrad*choose(cut, flow.grad(Functional(1), Identity(), ispgrad),
                         fhigh.grad(Functional(1), Identity(), ispgrad))(x);
  }
  void grad(const GridDescription &, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    const int N = data.cols()*data.rows();
    if (outpgrad) {
      for (int i=0; i<N; i++) {
        double x = data[i];
        if (x<cut) {
          (*outgrad)[i] += ingrad[i]*flow.derive(x);
          (*outpgrad)[i] += ingrad[i]*flow.derive(x);
        } else {
          (*outgrad)[i] += ingrad[i]*fhigh.derive(x);
           (*outpgrad)[i] += ingrad[i]*fhigh.derive(x);
        }
      }
    } else {
      for (int i=0; i<N; i++) {
        double x = data[i];
        (*outgrad)[i] += ingrad[i] * ((x<cut) ? flow.derive(x) : fhigh.derive(x));
      }
    }
  }

  Expression printme(const Expression &x) const {
    return funexpr("choose", cut, flow.printme(x), fhigh.printme(x));
  }
private:
  double cut;
  Functional flow, fhigh;
};

Functional choose(double cut, const Functional &lower, const Functional &higher) {
  return Functional(new Choose(cut, lower, higher));
}

static const double min_log_arg = 1e-90;
static const double slope = log(min_log_arg);
static const double min_e = min_log_arg*log(min_log_arg) - min_log_arg;

Functional find_nQ(double Tin) {
  Functional kT = Functional(Tin, "kT");
  Functional mass = 18*1822.8885; // FIXME: molecular weight of water
  // Note:  hbar is one in atomic units, yay!
  // nQ = (m*kT/(2 pi hbar^2))^3/2
  return PowAndHalf(3)(mass*kT/(2*M_PI));
}

Functional find_dnQ_dT(double Tin) {
  Functional kT = Functional(Tin, "kT");
  Functional mass = 18*1822.8885; // FIXME: molecular weight of water
  // Note:  hbar is one in atomic units, yay!
  // nQ = (m*kT/(2 pi hbar^2))^3/2
  return 1.5*PowAndHalf(3)(mass/(2*M_PI))*sqrt(kT);
}

Functional IdealGas(double Tin) {
  Functional n = Identity().set_name("n");
  Functional nQ = find_nQ(Tin);
  Functional T = Functional(Tin, "kT");
  return (T*choose(min_log_arg,  ((n-min_log_arg)*slope + min_e), (n*log(n/nQ)-n))).set_name("ideal gas");
}


Functional IdealGasOfVeff(double Tin) {
  Functional Veff = Identity().set_name("Veff");
  Functional kT = Functional(Tin, "kT");
  Functional n = exp(-Veff/kT);
  Functional nQ = find_nQ(Tin);
  return (-(Veff + kT*log(nQ) + kT)*n).set_name("ideal_gas");
}

Functional EntropyOfIdealGasOfVeff(double Tin) {
  Functional Veff = Identity().set_name("Veff");
  Functional kT = Functional(Tin, "kT");
  Functional n = exp(-Veff/kT);
  Functional nQ = find_nQ(Tin);
  Functional dnQ_dT = find_dnQ_dT(Tin);
  // The following is also known as the Sackur-Tetrode equation
  return (-(Veff/kT + log(nQ) + 2.5)*n).set_name("dideal_gas_dT");
}
