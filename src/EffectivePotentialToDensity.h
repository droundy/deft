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
    return (minusbeta*data).cwise().exp();
  }

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  void grad(const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += (minusbeta*ingrad).cwise()*((minusbeta*data).cwise().exp());
    if (outpgrad) {
      *outpgrad += (minusbeta*ingrad).cwise()*((minusbeta*data).cwise().exp());
    }
  }
private:
  double minusbeta; // -1/kT
};
