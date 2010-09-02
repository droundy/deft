// -*- mode: C++; -*-

#pragma once

#include "Functional.h"
#include <stdio.h>
#include <math.h>

class Minimizer;

Minimizer Downhill(Functional f, const GridDescription &gdin, VectorXd *data, double viscosity=0.1);
Minimizer PreconditionedDownhill(Functional f, const GridDescription &gdin, VectorXd *data, double viscosity=0.1);
Minimizer MaxIter(int maxiter, Minimizer);
Minimizer Precision(double err, Minimizer);

class MinimizerInterface {
public: // yuck, this shouldn't be public!
  Functional f;
  VectorXd *x; // Note that we don't own this data!
  int iter;
  GridDescription gd;
public:
  MinimizerInterface(Functional myf, const GridDescription &gdin, VectorXd *data)
    : f(myf), x(data), gd(gdin), last_grad(0), last_pgrad(0) {
    iter = 0;
  }
  virtual ~MinimizerInterface() {
    delete last_grad;
    delete last_pgrad;
  }
  virtual void minimize(Functional newf, const GridDescription &gdnew, VectorXd *newx = 0) {
    f = newf;
    iter = 0;
    invalidate_cache();
    if (newx) x = newx;
    gd = gdnew;
  }

  // improve_energy returns false if the energy is fully converged
  // (i.e. it didn't improve), and there is no reason to call this
  // minimizer any more.  Thus improve_energy can be naturally used as
  // the conditional in a while or for loop.
  virtual bool improve_energy(bool verbose = false) = 0;

  // The print_info function should be called at each iteration,
  // unless verbose is set to false in improve_energy.  But you can
  // also call it manually.
  virtual void print_info(const char *prefix = "") const;

  // energy returns the current energy.
  double energy() const { return f(gd, *x); }
  const VectorXd &grad() const {
    if (!last_grad) {
      last_grad = new VectorXd(*x); // hokey
      last_grad->setZero(); // Have to remember to zero it out first!
      f.grad(gd, *x, last_grad);
    }
    return *last_grad;
  }
  const VectorXd &pgrad() const {
    if (!last_pgrad) {
      last_pgrad = new VectorXd(*x); // hokey
      last_pgrad->setZero(); // Have to remember to zero it out first!
      if (!last_grad) last_grad = new VectorXd(*x); // hokey
      last_grad->setZero(); // Have to remember to zero it out first!
      f.grad(gd, *x, last_grad, last_pgrad);
    }
    return *last_pgrad;
  }

  // Note that we're changing the position x.
  void invalidate_cache() {
    delete last_grad;
    last_grad = 0;
    delete last_pgrad;
    last_pgrad = 0;
  }
private:
  mutable VectorXd *last_grad, *last_pgrad;
};

class Minimizer : public MinimizerInterface {
public:
  // Handle reference counting so we can pass these things around freely...
  explicit Minimizer(MinimizerInterface *p) // allocate a new counter
    : MinimizerInterface(p->f, p->gd, p->x), itsCounter(0) {
    itsCounter = new counter(p);
  }
  ~Minimizer() {
    release();
  }
  Minimizer(const Minimizer& r) : MinimizerInterface(r) {
    acquire(r.itsCounter);
  }
  Minimizer& operator=(const Minimizer& r) {
    if (this != &r) {
      release();
      acquire(r.itsCounter);
    }
    return *this;
  }
  void minimize(Functional newf, const GridDescription &gdnew, VectorXd *newx = 0) {
    MinimizerInterface::minimize(newf, gdnew, newx);
    itsCounter->ptr->minimize(newf, gdnew, newx);
  }
  bool improve_energy(bool verbose = false) {
    return itsCounter->ptr->improve_energy(verbose);
  }
  void print_info(const char *prefix = "") const {
    return itsCounter->ptr->print_info(prefix);
  }

private:
  struct counter {
    counter(MinimizerInterface* p) : ptr(p), count(1) {}
    MinimizerInterface* ptr;
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
      if (--itsCounter->count <= 0) {
        delete itsCounter->ptr;
        delete itsCounter;
      }
      itsCounter = 0;
    }
  }
};

// A base class handling a bit of busy work dealing with modifiers of
// minimizers.
class MinimizerModifier : public MinimizerInterface {
protected:
  Minimizer min;
public:
  MinimizerModifier(Minimizer m)
    : MinimizerInterface(m.f, m.gd, m.x), min(m) {}
  ~MinimizerModifier() { }
  void minimize(Functional newf, const GridDescription &gdnew, VectorXd *newx = 0) {
    MinimizerInterface::minimize(newf, gdnew, newx);
    min.minimize(newf, gdnew, newx);
  }

  bool improve_energy(bool verbose = false) {
    iter++;
    return min.improve_energy(verbose);
  }
  void print_info(const char *prefix="") const {
    return min.print_info(prefix);
  }
};

inline bool better(double a, double b) {
  return a < b || isnan(b);
}
