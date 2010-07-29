// -*- mode: C++; -*-

#pragma once

#include "FieldFunctional.h"

class EffectivePotentialToDensity : public FieldFunctional {
public:
  EffectivePotentialToDensity(double temp) {
    minusbeta = -1.0/temp;
  }

  // A functional mapping one field onto another...
  VectorXd operator()(const VectorXd &data) const {
    VectorXd out(data);
    for (int i=0; i<out.rows(); i++) {
      out[i] = exp(minusbeta*out[i]);
    }
    return out;
    //return (minusbeta*data).cwise().exp();
  }

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  void grad(const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    for (int i=0; i<data.rows(); i++) {
      (*outgrad)[i] += minusbeta*ingrad[i]*exp(minusbeta*data[i]);
    }
    if (outpgrad) *outpgrad += *outgrad;
    //*outgrad += (minusbeta*ingrad).cwise()*((minusbeta*data).cwise().exp());
    //if (outpgrad) {
    // *outpgrad += (minusbeta*ingrad).cwise()*((minusbeta*data).cwise().exp());
    //}
  }
private:
  double minusbeta; // -1/kT
};
