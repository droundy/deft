// -*- mode: C++; -*-

#pragma once

#include "Grid.h"

class FieldFunctionalInterface {
public:
  virtual ~FieldFunctionalInterface() {}

  // A functional mapping one field onto another...
  virtual VectorXd transform(const GridDescription &gd, const VectorXd &data) const = 0;
  virtual double transform(double) const = 0;

  // This computes the gradient of the functional, given a gradient of
  // its output field (i.e. it applies the chain rule).
  virtual void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
                    VectorXd *outgrad, VectorXd *outpgrad) const = 0;
};

class FieldFunctional {
public:
  // Handle reference counting so we can pass these things around freely...
  explicit FieldFunctional(FieldFunctionalInterface* p = 0, const char *name = 0) // allocate a new counter
    : itsCounter(0) {
    if (p) {
      itsCounter = new counter(p);
      itsCounter->name = name;
    }
    mynext = 0;
  }
  ~FieldFunctional() {
    release();
    delete mynext;
  }
  FieldFunctional(const FieldFunctional& r) {
    acquire(r.itsCounter);
    if (r.mynext) mynext = new FieldFunctional(*r.mynext);
    else mynext = 0;
  }
  FieldFunctional& operator=(const FieldFunctional& r) {
    if (this != &r) {
      release();
      acquire(r.itsCounter);
    }
    if (r.mynext) mynext = new FieldFunctional(*r.mynext);
    return *this;
  }

  FieldFunctional operator()(const FieldFunctional &) const;
  FieldFunctional operator+(const FieldFunctional &x) const {
    FieldFunctional out(*this);
    return out += x;
  }
  FieldFunctional operator+=(const FieldFunctional &x) {
    if (mynext) *mynext += x;
    else mynext = new FieldFunctional(x);
    return *this;
  }
  FieldFunctional operator-(const FieldFunctional &) const;
  FieldFunctional operator/(const FieldFunctional &) const;
  FieldFunctional operator*(const FieldFunctional &) const;
  FieldFunctional operator*(double) const;
  VectorXd justMe(const GridDescription &gd, const VectorXd &data) const {
    return itsCounter->ptr->transform(gd, data);
  }
  VectorXd operator()(const GridDescription &gd, const VectorXd &data) const {
    VectorXd out = itsCounter->ptr->transform(gd, data);
    if (mynext) out += (*mynext)(gd, data);
    return out;
  }
  VectorXd operator()(const Grid &g) const {
    VectorXd out = itsCounter->ptr->transform(g.description(), g);
    if (mynext) out += (*mynext)(g.description(), g);
    return out;
  }
  double operator()(double data) const {
    double out = itsCounter->ptr->transform(data);
    if (mynext) out += (*mynext)(data);
    return out;
  }
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    itsCounter->ptr->grad(gd, data, ingrad, outgrad, outpgrad);
    if (mynext) mynext->grad(gd, data, ingrad, outgrad, outpgrad);
  }
  const char *get_name() const { return itsCounter->name; }
  FieldFunctional set_name(const char *n) { itsCounter->name = n; return *this; }
  FieldFunctional *next() const {
    return mynext;
  }
  mutable double last_energy;
private:
  FieldFunctional *mynext;
  struct counter {
    counter(FieldFunctionalInterface* p = 0, unsigned c = 1) : ptr(p), count(c) {}
    FieldFunctionalInterface* ptr;
    unsigned count;
    const char *name;
    FieldFunctional *next;
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

inline FieldFunctional operator*(double x, const FieldFunctional &f) {
  return f*x;
}
FieldFunctional operator-(double x, const FieldFunctional &f);

FieldFunctional log(const FieldFunctional &);
FieldFunctional sqr(const FieldFunctional &);
