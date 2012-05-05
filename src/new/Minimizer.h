// -*- mode: C++; -*-

#pragma once

#include "Functional.h"
#include <stdio.h>
#include <math.h>

class Minimizer {
public:
  Minimizer(const Functional *myf, Vector *data)
    : f(myf), x(data), last_grad(0), last_pgrad(0) {
    iter = 0;
  }
  ~Minimizer() {
    delete last_grad;
    delete last_pgrad;
  }
  void minimize(const Functional *newf, Vector *newx = 0) {
    f = newf;
    iter = 0;
    invalidate_cache();
    if (newx) x = newx;
  }

  // improve_energy returns false if the energy is fully converged
  // (i.e. it didn't improve), and there is no reason to call this
  // minimizer any more.  Thus improve_energy can be naturally used as
  // the conditional in a while or for loop.
  bool improve_energy(bool verbose = false);

  // The print_info function should be called at each iteration,
  // unless verbose is set to false in improve_energy.  But you can
  // also call it manually.
  void print_info(const char *prefix = "") const;

  // energy returns the current energy.
  double energy() const {
    if (last_energy == 0) {
      last_energy = new double(f->energy(*x));
    }
    return *last_energy;
  }
  const Vector &grad() const {
    if (!last_grad) {
      last_grad = new Vector(f->grad(*x)); // a little sloppy here...
    }
    return *last_grad;
  }
  const Vector &pgrad() const {
    if (!last_pgrad) {
      if (f->have_preconditioner()) {
        invalidate_cache();
        EnergyGradAndPrecond foo = f->energy_grad_and_precond(*x, false, 0);
        last_energy = new double(foo.energy);
        last_grad = new Vector(foo.grad);
        last_pgrad = new Vector(foo.precond);
      } else {
        grad();
        last_pgrad = last_grad;
      }
    }
    return *last_pgrad;
  }

  // Note that we're changing the position x.
  void invalidate_cache() const {
    if (last_pgrad != last_grad) delete last_pgrad;
    last_pgrad = 0;
    delete last_grad;
    last_grad = 0;
    delete last_energy;
    last_energy = 0;
  }
private:
  const Functional *f;
  Vector *x; // Note that we don't own this data!
  int iter, maxiter;
  double precision;

  mutable Vector *last_grad, *last_pgrad;
  mutable double *last_energy;

  bool use_preconditioning;
  bool use_conjugate_gradient;
  bool check_conjugacy;

  double stepsize, orig_stepsize;
  Vector *direction, *oldgrad;

  double oldgradsqr;

};
