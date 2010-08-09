// -*- mode: C++; -*-

#pragma once

#include <eigen2/Eigen/Core>

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class FieldFunctionalInterface {
public:
  // A functional mapping one field onto another...
  virtual VectorXd operator()(const VectorXd &data) const = 0;
  virtual double operator()(double) const = 0;

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  virtual void grad(const VectorXd &data, const VectorXd &ingrad,
                    VectorXd *outgrad, VectorXd *outpgrad) const = 0;
};

class FieldFunctional : public FieldFunctionalInterface {
public:
  // Handle reference counting so we can pass these things around freely...
  explicit FieldFunctional(FieldFunctionalInterface* p = 0) // allocate a new counter
    : itsCounter(0) {
    if (p) itsCounter = new counter(p);
  }
  ~FieldFunctional() { release(); }
  FieldFunctional(const FieldFunctional& r) {
    acquire(r.itsCounter);
  }
  FieldFunctional& operator=(const FieldFunctional& r) {
    if (this != &r) {
      release();
      acquire(r.itsCounter);
    }
    return *this;
  }

  VectorXd operator()(const VectorXd &data) const {
    return (*itsCounter->ptr)(data);
  }
  double operator()(double data) const {
    return (*itsCounter->ptr)(data);
  }
  void grad(const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    itsCounter->ptr->grad(data, ingrad, outgrad, outpgrad);
  }
private:
  struct counter {
    counter(FieldFunctionalInterface* p = 0, unsigned c = 1) : ptr(p), count(c) {}
    FieldFunctionalInterface* ptr;
    unsigned count;
  };
  counter *itsCounter;
  void acquire(counter* c) {
    // increment the count
    itsCounter = c;
    if (c) (c->count)++;
  }
  void release() {
    // decrement the count, delete if it is 0
    if (itsCounter) {
      if (--itsCounter->count == 0) {
        delete itsCounter->ptr;
        delete itsCounter;
      }
      itsCounter = 0;
    }
  }
};
