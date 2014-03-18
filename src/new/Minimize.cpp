#include "Minimize.h"
#include <math.h>
#include <float.h>

inline bool better(double a, double b) {
  return a < b || b != b;
}

bool Minimize::improve_energy(Verbosity v) {
  iter++;
  if (iter >= maxiter) {
    if (v >= verbose) {
      printf("We reached the maximum number of iterations: %d with uncertainty %g remaining.\n",
             maxiter, error_estimate);
      fflush(stdout);
    }
    return false;
  }
  //printf("I am running ConjugateGradient::improve_energy\n");
  const double E0 = energy(v);
  const double old_deltaE = deltaE;
  if (E0 != E0) {
    // There is no point continuing, since we're starting with a NaN!
    // So we may as well quit here.
    if (v >= verbose) {
      printf("The initial energy is a NaN, so I'm quitting early.\n");
      fflush(stdout);
    }
    return false;
  }
  //f->run_finite_difference_test("functional");
  double gdotd;
  {
    const Vector pg = pgrad(v);
    const Vector g = grad();
    // Let's immediately free the cached gradient stored internally!
    invalidate_cache();

    // Note: my notation vaguely follows that of
    // [wikipedia](http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method).
    // I use the Polak-Ribiere method, with automatic direction reset.
    // Note that we could save some memory by using Fletcher-Reeves, and
    // it seems worth implementing that as an option for
    // memory-constrained problems (then we wouldn't need to store oldgrad).
    if (v >= min_details) {
      printf("\t\tnorm of gradient is %g\n", g.norm());
      //printf("\t\toldgrad size is %d\n", oldgrad.get_size());
      //printf("\t\tnorm(g - oldgrad) is %g\n", (g-oldgrad).norm());
    }
    double beta = pg.dot(g - oldgrad)/oldgradsqr;
    oldgrad = g;
    if (beta < 0 || beta != beta || oldgradsqr == 0) beta = 0;
    oldgradsqr = pg.dot(g);
    direction = -pg + beta*direction;
    gdotd = g.dot(direction);
    if (gdotd > 0) {
      direction = -g; // If our direction is uphill, reset to gradient.
      if (v >= verbose) printf("\t\treset to gradient, since g*d = %g < 0\n", gdotd);
      gdotd = g.dot(direction);
    }
    // g and pg will be destructed here.
  }

  {
    // Now we will do the line minimization... this is a bit
    // complicated.  We want to use as few steps as possible, but also
    // want to make sure we improve the energy at least a little bit.
    //if (v >= min_details) printf("\t\tInitial stepsize is %g\n", step);
  
    const double slope = -gdotd;
    if (v >= min_details) {
      printf("\t\tQuad: s0 = %25.15g  E0 = %25.15g", 0.0, E0);
      fflush(stdout);
    }
    if (E0 != E0) {
      // There is no point continuing, since we're starting with a NaN!
      // So we may as well quit here.
      if (v >= verbose) {
        printf(" which is a NaN, so I'm quitting early.\n");
        fflush(stdout);
      }
      return false;
    }
    if (E0 >= DBL_MAX || E0 <= -DBL_MAX) {
      // There is no point continuing, since we've got an infinite result.  :(
      // So we may as well quit here.
      if (v >= verbose) {
        printf(" which is infinite, so I'm quitting early.\n");
        fflush(stdout);
      }
      return false;
    }
    if (slope != slope) {
      // The slope here is a NaN, so there is no point continuing!
      // So we may as well quit here.
      if (v >= verbose) {
        printf(", but the slope is a NaN, so I'm quitting early.\n");
        fflush(stdout);
      }
      return false;
    }
    if (slope == 0) {
      // When the slope is precisely zero, there's no point continuing,
      // as the gradient doesn't have any information about how to
      // improve things.  Or possibly we're at the minimum to numerical
      // precision.
      if (v >= verbose) {
        printf("\nslope = %g which means we've arrived at the minimum!\n", slope);
        fflush(stdout);
      }
      return false;
    }

    if (slope*step < 0) {
      if (v >= verbose) printf("\n\t\tSwapping sign of step with slope %g...\n\t\t", slope);
      step *= -1;
    }
    // Here we're going to keep halving step1 until it's small enough
    // that the energy drops a tad.  This hopefully will give us a
    // reasonable starting guess.
    double step1 = step;
    *f += step1*direction; // Step away a bit...
    double Etried = -137;
    do {
      step1 *= 0.5;
      *f -= step1*direction; // and then step back a bit less...
      invalidate_cache();
      if (energy(v) == Etried) {
        if (v >= verbose) {
          printf("\tThis is silly in QuadraticLineMinimizeType::improve_energy: %g (%g vs %g)\n",
                 step1, energy(v), Etried);
          //Grid foo(gd, *x);
          //f.run_finite_difference_test("In QuadraticLineMinimizeType", kT, foo, &direction);
          //fflush(stdout);
        }
        break;
      }
      Etried = energy(v);
    } while (better(E0,energy(v)));
  
    const double E1 = energy(v);

    // Do a parabolic extrapolation with E0, E1, and gd to get curvature
    // and new step
    const double curvature = 2.0*(E1-E0+step1*slope)/(step1*step1);
    double step2 = slope/curvature;
    if (v >= min_details) {
      //printf("   E1 = %25.15g\n", E1);
      printf("\n\t\tQuad: s1 = %25.15g  E1 = %25.15g\n", step1, E1);
      //printf("\t\t                            Predicted E1 = %25.15g\n", E0 - step1*slope);
      printf("\t\tQuad: slope = %14.7g  curvature = %14.7g\n", slope, curvature);
      fflush(stdout);
    }
 
    if (curvature <= 0.0) {
      if (v >= verbose) {
        printf("\t\tCurvature has wrong sign... %g\n", curvature);
        fflush(stdout);
      }
      if (better(E0, E1)) {
        // It doesn't look to be working, so let's panic!
        invalidate_cache();
        *f -= step1*direction;
        if (v >= verbose) printf("\t\tQuadratic linmin not working properly!!!\n");
        //assert(!"Quadratic linmin not working well!!!\n");
        return false;
      } else {
        // It looks like we aren't moving far enough, so let's try
        // progressively doubling how far we're going.
        double Ebest = E0;
        while (energy(v) <= Ebest) {
          Ebest = energy(v);
          if (Ebest != E1 && v >= min_details)
            printf("\t\tQuad: si = %25.15g  Ei = %25.15g\n", step1, Ebest);
          *f += step1*direction;
          invalidate_cache();
          step1 *= 2;
        }
        if (v >= min_details) printf("\t\tQuad: sb = %25.15g  Eb = %25.15g\n", step1, energy(v));
        step1 *= 0.5;
        *f -= step1*direction;
        step = step1;
        invalidate_cache();
      }
    } else if (E1 == E0) {
      if (v >= verbose) {
        printf("\t\tNo change in energy... step1 is %g, direction has mag %g\n",
               step1, direction.norm());
        fflush(stdout);
      }
      do {
        invalidate_cache();
        *f -= step1*direction;
        step1 *= 2;
        *f += step1*direction;
      } while (energy(v) == E0);
      if (energy(v) > E0) {
        *f -= step1*direction;
        invalidate_cache();
        step1 = 0;
        printf("\t\tQuad: failed to find any improvement!  :(\n");
      } else if (energy(v) != energy(v)) {
        printf("\t\tQuad: found NaN with smallest possible step!!!\n");
        *f -= step1*direction;
        invalidate_cache();
        step1 = 0;
      }
    } else {
      step = step2; // output the stepsize for later reuse.
      invalidate_cache();
      *f += (step2-step1)*direction; // and move to the expected minimum.
      if (v >= min_details) {
        printf("\t\tQuad: s2 = %25.15g  E2 = %25.15g\n", step2, energy(v));
        fflush(stdout);
      }
      
      // Check that the energy did indeed drop!  FIXME: this may do
      // extra energy calculations, since it's not necessarily shared
      // with the driver routine!
      if (better(E1, energy(v)) && better(E1, E0)) {
        // The first try was better, so let's go with that one!
        if (v >= verbose) printf("\t\tGoing back to the first try...\n");
        invalidate_cache();
        *f -= (step2-step1)*direction;
        step = step1;
      } else if (energy(v) > E0) {
        const double E2 = energy(v);
        invalidate_cache();
        *f -= step2*direction;
        if (v >= verbose) {
          printf("\t\tQuadratic linmin not working well!!! (E2 = %14.7g, E0 = %14.7g)\n",
                 E2, energy(v));
          fflush(stdout);
        }
        return false;
        //assert(!"Quadratic linmin not working well!!!\n");
      }
    }
  }
  // Finished with the line minimization!

  if (v >= verbose) {
    pgrad(); // call pgrad here, since it might be more efficient to
             // compute everything at once, rather than first
             // computing just the grad...
    printf("\t\tfinal stepsize: = %g\n", step);
    if (v >= min_details) printf("\t\tgrad*direction = %g\n", grad().dot(direction)/gdotd);
    print_info("");
  }

  // At this point, we start work on estimating how close we are to
  // being adequately converged.
  const double newE = energy(v);
  deltaE = newE - E0;
  
  const double w = 0.01; // weighting for windowed average
  dEdn = (1-w)*dEdn + w*deltaE;
  dEdn = max(deltaE, old_deltaE);

  const double new_log_dEdn_ratio_average =
    (deltaE && old_deltaE) ? log(fabs(deltaE/old_deltaE)) : 0;
  log_dEdn_ratio_average = (1-w)*log_dEdn_ratio_average + w*new_log_dEdn_ratio_average;

  const double dEdn_ratio_average = exp(log_dEdn_ratio_average);
  // We assume below an exponential again...
  double error_guess = fabs(max(dEdn,deltaE)/log_dEdn_ratio_average);
  if (dEdn_ratio_average >= 1) {
    // We aren't converging at all! We'll just fudge a guess here,
    // adding on a bit of precision to make sure we don't stop early.
    error_guess = precision + fabs(newE);
  }

  error_estimate = 2*error_guess; // Just a bit of paranoia...

  if (known_true_energy && v >= verbose) {
    printf("True error is:  %10g\n", energy() - known_true_energy);
  }
  if (deltaE == 0 && old_deltaE == 0) {
    if (v >= verbose) printf("We got no change twice in a row, so we're done!\n");
    if (iter < miniter) return true;
    return false;
  }
  if (!(error_estimate < precision) && !(error_estimate < relative_precision*fabs(newE))) {
    if (dEdn_ratio_average < 1) {
      const double itersleft = log(precision/error_guess)/log(dEdn_ratio_average);
      if (v >= verbose) {
        if (precision > 0) {
          printf("Error estimate is %10g ... %.1f iterations remaining.\n",
                 error_estimate, itersleft);
        } else {
          printf("Error estimate is %10g.\n", error_estimate);
        }
        printf("\n");
      }
    }
    return true;
  } else {
    if (v >= verbose) printf("Converged with precision of %g!\n", error_estimate);
    if (iter < miniter) return true;
    return false;
  }
}

void Minimize::print_info(const char *prefix, bool with_iteration) const {
  if (with_iteration) printf("\n%s==== Iteration %d ====\n", prefix, iter);
  f->printme(prefix);
  printf("%s total energy =", prefix);
  print_double("", energy(silent));
  printf("\n\n");
}
