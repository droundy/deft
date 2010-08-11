// -*- mode: C++; -*-

#pragma once

#include "Functional.h"

class Minimizer {
protected:
  Functional f;
  VectorXd *x; // Note that we don't own this data!
  int iter;
public:
  Minimizer(Functional myf, VectorXd *data)
    : f(myf), x(data), last_grad(0), last_pgrad(0) {
    iter = 0;
  }
  ~Minimizer() {
    delete last_grad;
    delete last_pgrad;
  }
  virtual void minimize(Functional newf, VectorXd *newx = 0) {
    f = newf;
    iter = 0;
    invalidate_cache();
    if (newx) x = newx;
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
  double energy() const { return f(*x); }
  const VectorXd &grad() const {
    if (!last_grad) {
      last_grad = new VectorXd(*x); // hokey
      last_grad->setZero(); // Have to remember to zero it out first!
      f.grad(*x, last_grad);
    }
    return *last_grad;
  }
  const VectorXd &pgrad() const {
    if (!last_pgrad) {
      last_pgrad = new VectorXd(*x); // hokey
      last_pgrad->setZero(); // Have to remember to zero it out first!
      if (!last_grad) {
        last_grad = new VectorXd(*x); // hokey
        last_grad->setZero(); // Have to remember to zero it out first!
      }
      f.grad(*x, last_grad, last_pgrad);
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
