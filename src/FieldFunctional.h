// -*- mode: C++; -*-

#pragma once

#include "ReciprocalGrid.h"

class FieldFunctional;

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
  virtual double grad(double data) const = 0;
  //virtual FieldFunctional grad(const FieldFunctional &ingrad) const = 0;

  virtual void print_summary(const char *prefix, double energy, const char *name) const;
};

template<typename Derived>
class ConvolveWith : public FieldFunctionalInterface {
public:
  ConvolveWith(const Eigen::MatrixBase<Derived> &x) : c(x) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    Grid out(gd, data);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= c;
    return recip.ifft();
  }
  double transform(double n) const {
    return n*c[0];
  }
  double grad(double) const {
    return c[0];
  }
  void grad(const GridDescription &gd, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= c;
    out = recip.ifft();
    *outgrad += out;
    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
private:
  Derived c;
};

class FieldFunctional {
public:
  // Handle reference counting so we can pass these things around freely...
  FieldFunctional(double); // This handles constants!
  explicit FieldFunctional(const VectorXd &); // This handles constant fields!
  template<typename Derived>
  explicit FieldFunctional(const MatrixBase<Derived> &x, const char *name) : itsCounter(0) {
    // This handles constant ephemeral fields!
    init(new ConvolveWith<Derived>(x), name);
  }
  explicit FieldFunctional(FieldFunctionalInterface* p = 0, const char *name = 0) // allocate a new counter
    : itsCounter(0) {
    init(p, name);
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
  double integral(const Grid &g) const {
    return integral(g.description(), g);
  }
  double integral(const GridDescription &gd, const VectorXd &data) const {
    // This does some extra work to save the energies of each term in
    // the sum.
    VectorXd fdata(justMe(gd, data));
    double e = gd.dvolume*fdata.sum();
    set_last_energy(e);
    FieldFunctional *nxt = next();
    while (nxt) {
      fdata += nxt->justMe(gd, data);
      double etot = gd.dvolume*fdata.sum();
      nxt->set_last_energy(etot - e);
      e = etot;
      nxt = nxt->next();
    }
    return e;
  }
  void integralgrad(const Grid &g, VectorXd *gr, VectorXd *pg=0) const {
    grad(g.description(), g, g.description().dvolume*VectorXd::Ones(g.description().NxNyNz), gr, pg);
  }
  void integralgrad(const GridDescription &gd, const VectorXd &x, VectorXd *g, VectorXd *pg=0) const {
    grad(gd, x, gd.dvolume*VectorXd::Ones(gd.NxNyNz), g, pg);
  }
  double operator()(double data) const {
    double out = itsCounter->ptr->transform(data);
    if (mynext) out += (*mynext)(data);
    return out;
  }
  //FieldFunctional grad(const FieldFunctional &ingrad) const {
  //  return itsCounter->ptr->grad(ingrad);
  //}
  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    itsCounter->ptr->grad(gd, data, ingrad, outgrad, outpgrad);
    if (mynext) mynext->grad(gd, data, ingrad, outgrad, outpgrad);
  }
  double grad(double data) const {
    double out = itsCounter->ptr->grad(data);
    if (mynext) out += mynext->grad(data);
    return out;
  }
  const char *get_name() const { return itsCounter->name; }
  FieldFunctional set_name(const char *n) { itsCounter->name = n; return *this; }
  FieldFunctional *next() const {
    return mynext;
  }
  FieldFunctional set_last_energy(double e) const { itsCounter->last_energy = e; return *this; }

  void print_summary(const char *prefix, double energy, const char *name=0) const;
  // The following utility methods do not need to be overridden.  Its
  // return value is the total energy that is printed.
  double print_iteration(const char *prefix, int iter) const;
  // run_finite_difference_test returns false when the test fails.
  bool run_finite_difference_test(const char *testname,
                                  const Grid &data,
                                  const VectorXd *direction = 0) const;
private:
  void init(FieldFunctionalInterface *p, const char *name) {
    if (p) {
      itsCounter = new counter(p);
      itsCounter->name = name;
    }
    mynext = 0;
  }
  FieldFunctional *mynext;
  struct counter {
    counter(FieldFunctionalInterface* p = 0, unsigned c = 1) : ptr(p), count(c) {}
    FieldFunctionalInterface* ptr;
    unsigned count;
    const char *name;
    FieldFunctional *next;
    double last_energy;
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
FieldFunctional constrain(const Grid &, FieldFunctional);

// choose combines two local functionals, with which one is applied
// depending on the local value of the field.
FieldFunctional choose(double, const FieldFunctional &lower, const FieldFunctional &higher);
