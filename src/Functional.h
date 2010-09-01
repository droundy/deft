// -*- mode: C++; -*-

#pragma once

#include "FieldFunctional.h"
#include "Grid.h"
#include <stdint.h>
#include <stdio.h>

class FunctionalComposition;

class FunctionalInterface {
public:
  FunctionalInterface() {}
  virtual ~FunctionalInterface() {}

  // To implement a functional, you need to provide both an energy
  // method and a gradient method.

  // The gradient method, however, outputs its output in-place in two
  // VectorXd pointers.  If the second pointer is nonzero, you need to
  // also output a preconditioned gradient.  The two gradients should
  // be accumulated, so that the resulting VectorXd will be the sum of
  // the gradient plus any additional value.  This enables the easy
  // combination of different functionals.
  virtual void grad(const GridDescription &gd, const VectorXd &data,
                    VectorXd *, VectorXd *pgrad = 0) const = 0;

  // You may optionally define a print_summary method, which would
  // print something interesting to the screen.
  virtual void print_summary(const char *prefix, double energy) const;

  virtual double energy(const GridDescription &gd, const VectorXd &data) const = 0;
  virtual double energy(double data) const = 0;
};

class Functional {
public:
  // Handle reference counting so we can pass these things around freely...
  explicit Functional(FunctionalInterface* p = 0, const char *n = 0) // allocate a new counter
    : itsCounter(0) {
    if (p) {
      itsCounter = new counter(p);
      itsCounter->checksum = -1;
      itsCounter->name = n;
    }
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

  // The energy method is operator()
  double operator()(const GridDescription &gd, const VectorXd &data) const {
    uint32_t sum = compute_checksum(data);
    // FIXME: should also check that gd hasn't changed!
    if (sum != itsCounter->checksum || int(itsCounter->checksum) == -1) {
      itsCounter->last_energy = energy(gd, data);
      itsCounter->checksum = sum;
    }
    return itsCounter->last_energy;
  }
  double operator()(const Grid &g) const {
    uint32_t sum = compute_checksum(g);
    // FIXME: should also check that g.description() hasn't changed!
    if (sum != itsCounter->checksum || int(itsCounter->checksum) == -1) {
      itsCounter->last_energy = energy(g.description(), g);
      itsCounter->checksum = sum;
    }
    return itsCounter->last_energy;
  }
  // The following method is for a homogeneous density.
  double operator()(double data) const { return energy(data); }
  Functional operator()(const FieldFunctional &) const;

  double energy(double data) const {
    return itsCounter->ptr->energy(data);
  }
  double energy(const GridDescription &gd, const VectorXd &data) const {
    return itsCounter->ptr->energy(gd, data);
  }
  void grad(const GridDescription &gd, const VectorXd &data,
            VectorXd *g, VectorXd *pg = 0) const {
    itsCounter->ptr->grad(gd, data, g, pg);
  }
  void grad(const Grid &g, VectorXd *myg, VectorXd *mypg = 0) const {
    grad(g.description(), g, myg, mypg); // Just pull out the description...
  }

  // The following utility methods do not need to be overridden.
  void print_iteration(const char *prefix, int iter) const;
  // run_finite_difference_test returns false when the test fails.
  bool run_finite_difference_test(const char *testname,
                                  const Grid &data,
                                  const VectorXd *direction = 0) const;

  void print_summary(const char *prefix) const;

  void set_name(const char *n) { itsCounter->name = n; }
  const char *get_name() const { return itsCounter->name; }
private:
  struct counter {
    counter(FunctionalInterface* p = 0, unsigned c = 1, const char *n = 0) : ptr(p), count(c), name(n) {}
    FunctionalInterface* ptr;
    mutable uint32_t checksum;
    mutable double last_energy;
    unsigned count;
    const char *name;
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
  uint32_t compute_checksum(const VectorXd &m) const {
    uint32_t *mints = (uint32_t *)m.data();
    const int N = m.cols()*m.rows()*2;
    uint32_t sum = 0, othersum = -1;
    for (int i=0; i<N; i++) {
      sum ^= mints[i]+i; // FIXME: this is a really stupid checksum...
      othersum += (mints[i] << (i&7)) + (mints[i] >> (32-(i&7)));
    }
    return sum ^ othersum;
  }
};

Functional operator+(const Functional &, const Functional &);
Functional operator*(const Functional &, const Functional &);
Functional operator*(double, const Functional &);
inline Functional operator-(const Functional &f) {
  return (-1)*f;
}
extern Functional integrate;
Functional constrain(const Grid &, Functional);
