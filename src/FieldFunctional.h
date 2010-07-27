// -*- mode: C++; -*-

#pragma once

#include <eigen2/Eigen/Core>

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class FieldFunctional {
public:
  // A functional mapping one field onto another...
  virtual VectorXd operator()(const VectorXd &data) const = 0;

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  virtual void grad(const VectorXd &data, const VectorXd &ingrad,
                    VectorXd *outgrad, VectorXd *outpgrad) const = 0;
};
