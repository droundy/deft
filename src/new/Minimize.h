// -*- mode: C++; -*-

#pragma once

#include "new/NewFunctional.h"
#include "handymath.h"
#include <stdio.h>
#include <math.h>
#include <time.h>

const Verbosity min_details = chatty;

class Minimize {
public:
  Minimize(NewFunctional *myf) : f(myf) {
    iter = 0;
    maxiter = 10000000;
    miniter = 6;
    num_energy_calcs = 0;
    num_grad_calcs = 0;
    precision = 0;
    
    last_energy = 0;
    known_true_energy = 0;

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
    error_estimate = 0;
  }
  ~Minimize() {
    invalidate_cache();
  }
  void minimize(NewFunctional *newf) {
    f = newf;
    iter = 0;
    num_energy_calcs = 0;
    num_grad_calcs = 0;
    invalidate_cache();
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
  void set_miniter(int mn) {
    miniter = mn;
  }
  int get_iteration_count() const {
    return iter;
  }
  void set_known_true_energy(double e) {
    known_true_energy = e;
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
  void print_info(const char *prefix = "", bool with_iteration = true) const;

  // energy returns the current energy.
  double energy(Verbosity v = quiet) const {
    if (last_energy == 0) {
      const clock_t start = clock();

      last_energy = new double(f->energy());
      if (v >= louder(min_details)) { // we need to get really paranoid before we print each energy...
        const clock_t end = clock();
        if (end > start + 10) {
          printf("\n\t\t\t=== Energy calculation #%d (took %g seconds) ===\n",
                 num_energy_calcs, (end - double(start))/CLOCKS_PER_SEC);
        } else {
          // no point printing the time, since it probably isn't
          // accurate anyhow...
          printf("\n\t\t\t=== Energy calculation #%d ===\n", num_energy_calcs);
        }
        print_info("\t\t\t", false);
      }
      num_energy_calcs++;
    }
    return *last_energy;
  }
  const Vector &grad() const {
    if (!last_grad.get_size()) {
      last_grad = f->grad();
      num_grad_calcs++;
    }
    return last_grad;
  }
  const Vector &pgrad(Verbosity v = quiet) const {
    if (!last_pgrad.get_size()) {
      if (use_preconditioning && f->have_preconditioner()) {
        invalidate_cache();
        const clock_t start = clock();
        EnergyGradAndPrecond foo = f->energy_grad_and_precond();
        if (v >= louder(min_details)) { // we need to get really paranoid before we print each energy...
          const clock_t end = clock();
          if (end > start + 10) {
            printf("\n\t\t\t=== Energy and grad calculation #%d and #%d resp. (took %g seconds) ===\n",
                   num_energy_calcs, num_grad_calcs, (end - double(start))/CLOCKS_PER_SEC);
          } else {
            // no point printing the time, since it probably isn't
            // accurate anyhow...
            printf("\n\t\t\t=== Energy calculation #%d (grad #%d) ===\n", num_energy_calcs, num_grad_calcs);
          }
          print_info("\t\t\t", false);
        }
        num_energy_calcs++;
        num_grad_calcs++;
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
  NewFunctional *f;
  int iter, maxiter, miniter;

  mutable int num_energy_calcs, num_grad_calcs;
  mutable Vector last_grad, last_pgrad;
  mutable double *last_energy;

  bool use_preconditioning;
  bool use_conjugate_gradient;
  bool do_check_conjugacy;

  double step, orig_stepsize;
  Vector direction, oldgrad;

  double oldgradsqr;

  double precision, relative_precision, deltaE, dEdn, log_dEdn_ratio_average;
  double error_estimate;
  double known_true_energy; // used for checking how well the minimization is working
};
