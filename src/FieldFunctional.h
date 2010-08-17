// -*- mode: C++; -*-

#pragma once

#include "Grid.h"

class FieldFunctionalInterface {
public:
  virtual ~FieldFunctionalInterface() {}

  // A functional mapping one field onto another...
  virtual VectorXd operator()(const GridDescription &gd, const VectorXd &data) const = 0;
  virtual double operator()(double) const = 0;

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  virtual void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
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

  VectorXd operator()(const GridDescription &gd, const VectorXd &data) const {
    return (*itsCounter->ptr)(gd, data);
  }
  VectorXd operator()(const Grid &g) const {
    return (*this)(g.description(), g);
  }
  double operator()(double data) const {
    return (*itsCounter->ptr)(data);
  }
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    itsCounter->ptr->grad(gd, data, ingrad, outgrad, outpgrad);
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
