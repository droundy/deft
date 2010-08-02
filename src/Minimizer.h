// -*- mode: C++; -*-

#pragma once

#include "Functional.h"
#include "counted_ptr.h"

class Minimizer {
protected:
  counted_ptr<Functional> f;
  VectorXd *x; // Note that we don't own this data!
  int iter;
public:
  Minimizer(counted_ptr<Functional> myf, VectorXd *data)
    : f(myf), x(data), last_grad(*data) {
    iter = 0;
    energy_valid = grad_valid = false;
  }

  // improve_energy returns true if the energy is fully converged, and
  // there is no reason to call this minimizer any more.
  virtual bool improve_energy(bool verbose = false) = 0;

  // The print_info function should be called at each iteration,
  // unless verbose is set to false in improve_energy.  But you can
  // also call it manually.
  virtual void print_info(int iter) const;

  // energy returns the current energy.
  double energy() const {
    if (!energy_valid) {
      last_energy = (*f)(*x);
      energy_valid = true;
    }
    return last_energy;
  }
  const VectorXd &grad() const {
    if (!grad_valid) {
      last_grad.setZero(); // Have to remember to zero it out first!
      f->grad(*x, &last_grad);
      grad_valid = true;
    }
    return last_grad;
  }

  // Note that we're changing the position x.
  void invalidate_cache() {
    energy_valid = false;
    grad_valid = false;
  }
private:
  mutable double last_energy;
  mutable VectorXd last_grad;
  mutable bool energy_valid;
  mutable bool grad_valid;
};
