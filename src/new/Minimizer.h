// -*- mode: C++; -*-

#pragma once

#include "Functional.h"
#include <stdio.h>
#include <math.h>

class Minimizer {
public:
  Minimizer(const Functional *myf, Vector *data)
    : f(myf), x(data) {
    iter = 0;
    maxiter = 10000000;
    precision = 0;
    
    last_energy = 0;

    use_preconditioning = false;
    use_conjugate_gradient = false;
    do_check_conjugacy = false;

    step = 0.1;
    orig_stepsize = 0.1;

    oldgradsqr = 1.0;

    precision = 0;
    relative_precision = 1e-15;
    deltaE = 0;
    dEdn = 0;
    log_dEdn_ratio_average = 0;
  }
  ~Minimizer() {
    invalidate_cache();
  }
  void minimize(const Functional *newf, Vector *newx = 0) {
    f = newf;
    iter = 0;
    invalidate_cache();
    if (newx) x = newx;
  }

  // The following allow you to configure the algorithm used by the
  // minimizer.
  void set_precision(double p) {
    precision = p;
  }
  void set_relative_precision(double p) {
    relative_precision = p;
  }
  void set_maxiter(int mx) {
    maxiter = mx;
  }
  void precondition(bool u) {
    use_preconditioning = u;
  }
  void conjugate_gradient(bool u) {
    use_conjugate_gradient = u;
  }
  void check_conjugacy(bool u) {
    do_check_conjugacy = u;
  }

  // improve_energy returns false if the energy is fully converged
  // (i.e. it didn't improve), and there is no reason to call this
  // minimizer any more.  Thus improve_energy can be naturally used as
  // the conditional in a while or for loop.
  bool improve_energy(Verbosity verbose = quiet);

  // The print_info function should be called at each iteration,
  // unless verbose is set to quiet in improve_energy.  But you can
  // also call it manually.
  void print_info(const char *prefix = "") const;

  // energy returns the current energy.
  double energy(Verbosity v = quiet) const {
    if (last_energy == 0) {
      last_energy = new double(f->energy(*x, v));
    }
    return *last_energy;
  }
  const Vector &grad() const {
    if (!last_grad.get_size()) {
      last_grad = f->grad(*x);
    }
    return last_grad;
  }
  const Vector &pgrad() const {
    if (!last_pgrad.get_size()) {
      if (f->have_preconditioner()) {
        invalidate_cache();
        EnergyGradAndPrecond foo = f->energy_grad_and_precond(*x, false, 0);
        last_energy = new double(foo.energy);
        last_grad = foo.grad;
        last_pgrad = foo.precond;
      } else {
        grad();
        last_pgrad = last_grad;
      }
    }
    return last_pgrad;
  }

  // Note that we're changing the position x.
  void invalidate_cache() const {
    last_grad.free();
    last_pgrad.free();
    delete last_energy;
    last_energy = 0;
  }
private:
  const Functional *f;
  Vector *x; // Note that we don't own this data!
  int iter, maxiter;

  mutable Vector last_grad, last_pgrad;
  mutable double *last_energy;

  bool use_preconditioning;
  bool use_conjugate_gradient;
  bool do_check_conjugacy;

  double step, orig_stepsize;
  Vector direction, oldgrad;

  double oldgradsqr;

  double precision, relative_precision, deltaE, dEdn, log_dEdn_ratio_average;
};
