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

class EffectivePotentialToDensityType : public FieldFunctionalInterface {
public:
  EffectivePotentialToDensityType(double temp) {
    minusbeta = -1.0/temp;
  }

  // A functional mapping one field onto another...
  VectorXd operator()(const GridDescription &, const VectorXd &data) const {
    VectorXd out(data);
    for (int i=0; i<out.rows(); i++) {
      const double n = exp(minusbeta*out[i]);
      if (n != n) {
        printf("got NaN from out[%d] == %g\n", i, out[i]);
        assert(n == n);
      }
      out[i] = n;
    }
    return out;
    //return (minusbeta*data).cwise().exp();
  }
  double operator()(double Veff) const {
    return exp(minusbeta*Veff);
  }

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  void grad(const GridDescription &, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    for (int i=0; i<data.rows(); i++)
      (*outgrad)[i] += minusbeta*ingrad[i]*exp(minusbeta*data[i]);
    if (outpgrad) // precondition by not dividing by n! (big change)
      for (int i=0; i<data.rows(); i++) (*outpgrad)[i] += minusbeta*ingrad[i];
    //*outgrad += (minusbeta*ingrad).cwise()*((minusbeta*data).cwise().exp());
    //if (outpgrad) {
    // *outpgrad += (minusbeta*ingrad).cwise()*((minusbeta*data).cwise().exp());
    //}
  }
private:
  double minusbeta; // -1/kT
};

FieldFunctional EffectivePotentialToDensity(double temp) {
  return FieldFunctional(new EffectivePotentialToDensityType(temp));
}
