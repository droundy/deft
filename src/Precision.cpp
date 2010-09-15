// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Deft Authors
//
// Deft is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with deft; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include "Minimizer.h"
#include <stdio.h>

static inline double max(double a, double b) { return (a > b) ? a : b; }

// Impose a maximum number of iterations...
class PrecisionMinimizer : public MinimizerModifier {
protected:
  // dEdn is just the averaged slope, which may not be very
  // accurate if we're converging rapidly.
  double dEdn;
  // log_dEdn_ratio_average should be -gamma if
  //
  //     E(n) = Err(0) exp(-gamma*n) + E(infty)
  double log_dEdn_ratio_average;
  double error_estimate;
  double w;
  // The incremental change on the last iteration.
  double deltaE;

  double convergence_criterion;
  double itersleft;
public:
  PrecisionMinimizer(Minimizer m, double precision)
    : MinimizerModifier(m), convergence_criterion(precision) {
    error_estimate = max(1.0, 2*precision);
    dEdn = 0.0;
    log_dEdn_ratio_average = 0;
    deltaE = 0.0;
    w = 0.1;
  }
  ~PrecisionMinimizer() {}
  void minimize(FieldFunctional newf, const GridDescription &gdnew, VectorXd *newx = 0) {
    MinimizerModifier::minimize(newf, gdnew, newx);
  }

  bool improve_energy(bool verbose = false) {
    const double old_deltaE = deltaE;
    const double old_energy = MinimizerModifier::energy();
    bool not_done = MinimizerModifier::improve_energy(verbose);
    const double new_energy = MinimizerModifier::energy();
    deltaE = old_energy - new_energy;

    //const double old_dEdn = dEdn;
    dEdn = (1-w)*dEdn + w*deltaE;
    dEdn = max(deltaE, old_deltaE);

    //const double old_log_dEdn_ratio_average = log_dEdn_ratio_average;
    const double new_log_dEdn_ratio_average = (deltaE && old_deltaE) ? log(fabs(deltaE/old_deltaE)) : 0;
    log_dEdn_ratio_average = (1-w)*log_dEdn_ratio_average + w*new_log_dEdn_ratio_average;

    const double dEdn_ratio_average = exp(log_dEdn_ratio_average);
    // We assume below an exponential again...
    double error_guess = fabs(max(dEdn,deltaE)/log_dEdn_ratio_average);
    if (dEdn_ratio_average >= 1) {
      // We aren't converging at all! We'll just fudge a guess here,
      // adding on a bit of convergence_criterion to make sure we
      // don't stop early.
      error_guess = convergence_criterion + fabs(new_energy);
    }

    error_estimate = 2*error_guess; // Just a bit of paranoia...

    if (verbose) {
      //printf("dEdn       average is %10g ...\n", dEdn);
      //printf("dEdn_ratio         is %10g ...\n", deltaE/old_deltaE);
      //printf("dEdn_ratio_average is %10g ...\n", dEdn_ratio_average);
      //printf("error_guess        is %10g ...\n", error_guess);
    }
    if (deltaE == 0 && old_deltaE == 0) {
      if (verbose) printf("We got no change twice in a row, so we're done!\n");
      return false;
    }

    /*
    if (old_deltaE) {
      double old_log_rate_average = log(rate_average);
      double latest_log_rate_average = log(fabs(old_deltaE/deltaE));
      if (deltaE == 0.0) latest_log_rate_average = 0.0;
      double log_rate_average = (1-window)*old_log_rate_average +
        window*latest_log_rate_average;
      rate_average = exp(log_rate_average);
      
      const double err_window = 1.0 - 0.5/rate_average;
      if (deltaE > error_average)
        error_average = fabs(deltaE)/fabs(rate_average - 1.0);
      else if (deltaE > old_deltaE)
        error_average -= old_deltaE;
      else
        error_average -= deltaE;
      error_average *= 1.0 - err_window;
      error_average += err_window*fabs(deltaE)/fabs(rate_average - 1.0);
      if (deltaE == 0.0) deltaE = old_deltaE; // don't give up!
      itersleft = log(max(old_error_average, error_average)/
                      convergence_criterion)/log(rate_average);
      if (verbose) {
        printf("Error old del is %10g ...\n", old_deltaE);
        printf("Error delta E is %10g ...\n", deltaE);
        printf("Error rate    is %10g ...\n", rate_average);
        printf("Energy prediction %.14g ...\n", new_energy - error_average);
        if (convergence_criterion) {
          if (itersleft > 0)
            printf("Error estimate is %10g ... %.1f iterations remaining.\n",
                    error_average, itersleft);
          else printf("Error estimate is %10g ...\n",
                       max(error_average, old_error_average));
        } else {
          printf("Error estimate is %10g\n", error_average);
        }
      }
    } else if (deltaE == 0.0) {
      // we've got zero changes two times in a row, so we're done!!!
      error_average = 0.0;
    } else {
      // This is the first step... take a wild guess!
      error_average = fabs(deltaE) + convergence_criterion;
    }
    */
    if (not_done && !(error_estimate < convergence_criterion)) {
      if (dEdn_ratio_average < 1) {
        itersleft = log(convergence_criterion/error_guess)/log(dEdn_ratio_average);
        if (verbose) {
          if (convergence_criterion > 0) {
            printf("Error estimate is %10g ... %.1f iterations remaining.\n", error_estimate, itersleft);
          } else {
            printf("Error estimate is %10g.\n", error_estimate);
          }
        }
      }
      return true;
    } else {
      if (verbose) printf("Converged with precision of %g!\n", error_estimate);
      return false;
    }
  }
};

Minimizer Precision(double precision, Minimizer m) {
  return Minimizer(new PrecisionMinimizer(m, precision));
}
