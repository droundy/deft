// -*- mode: C++; -*-

#pragma once

#include "FieldFunctional.h"
#include <stdint.h>

class FunctionalComposition;

class FunctionalInterface {
public:
  FunctionalInterface() {
    checksum = -1;
  }

  // To implement a functional, you need to provide both an energy
  // method and a gradient method.

  // The energy method is operator()
  double operator()(const VectorXd &data) const {
    uint32_t sum = compute_checksum(data);
    if (sum != checksum || int(checksum) == -1) {
      last_energy = energy(data);
      checksum = sum;
    }
    return last_energy;
  }
  // The following method is for a homogeneous density.
  double operator()(double data) const { return energy(data); }

  // The gradient method, however, outputs its output in-place in two
  // VectorXd pointers.  If the second pointer is nonzero, you need to
  // also output a preconditioned gradient.  The two gradients should
  // be accumulated, so that the resulting VectorXd will be the sum of
  // the gradient plus any additional value.  This enables the easy
  // combination of different functionals.
  virtual void grad(const VectorXd &data,
                    VectorXd *, VectorXd *pgrad = 0) const = 0;

  // You may optionally define a print_summary method, which would
  // print something interesting to the screen.
  virtual void print_summary(const char *prefix, const VectorXd &data) const;

  // The following utility methods do not need to be overridden.
  void print_iteration(const char *prefix, const VectorXd &data, int iter) const;
  // run_finite_difference_test returns false when the test fails.
  bool run_finite_difference_test(const char *testname,
                                  const VectorXd &data,
                                  const VectorXd *direction = 0) const;
protected:
  mutable double last_energy;
  virtual double energy(const VectorXd &data) const = 0;
  virtual double energy(double data) const = 0;
private:
  mutable uint32_t checksum;
  uint32_t compute_checksum(const VectorXd &m) const {
    uint32_t *mints = (uint32_t *)m.data();
    const int N = m.cols()*m.rows()*2;
    uint32_t sum = 0;
    for (int i=0; i<N; i++) sum ^= mints[i]+i; // FIXME: this is a really stupid checksum...
    return sum;
  }
};

class Functional : public FunctionalInterface {
public:
  // Handle reference counting so we can pass these things around freely...
  explicit Functional(FunctionalInterface* p = 0) // allocate a new counter
    : itsCounter(0) {
    if (p) itsCounter = new counter(p);
  }
  ~Functional() { release(); }
  Functional(const Functional& r) {
    acquire(r.itsCounter);
  }
  Functional& operator=(const Functional& r) {
    if (this != &r) {
      release();
      acquire(r.itsCounter);
    }
    return *this;
  }

  double energy(double data) const {
    return (*itsCounter->ptr)(data);
  }
  double energy(const VectorXd &data) const {
    return (*itsCounter->ptr)(data);
  }
  void grad(const VectorXd &data,
            VectorXd *g, VectorXd *pg = 0) const {
    itsCounter->ptr->grad(data, g, pg);
  }
  void print_summary(const char *prefix, const VectorXd &data) const {
    itsCounter->ptr->print_summary(prefix, data);
  }
private:
  struct counter {
    counter(FunctionalInterface* p = 0, unsigned c = 1) : ptr(p), count(c) {}
    FunctionalInterface* ptr;
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

Functional operator+(const Functional &, const Functional &);
Functional compose(const Functional &a, const FieldFunctional &b);
