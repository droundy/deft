#include "Minimizer.h"
#include <math.h>

inline bool better(double a, double b) {
  return a < b || isnan(b);
}

Verbosity quieter(Verbosity v) {
  return Verbosity(v-1);
}

Verbosity quietest(Verbosity v) {
  return Verbosity(v-2);
}

const Verbosity min_details = chatty;

bool Minimizer::improve_energy(Verbosity v) {
  iter++;
  //printf("I am running ConjugateGradient::improve_energy\n");
  const double E0 = energy(quietest(v));
  if (isnan(E0)) {
    // There is no point continuing, since we're starting with a NaN!
    // So we may as well quit here.
    if (v >= verbose) {
      printf("The initial energy is a NaN, so I'm quitting early.\n");
      fflush(stdout);
    }
    return false;
  }
  double gdotd;
  {
    const Vector g = -grad();
    // Let's immediately free the cached gradient stored internally!
    invalidate_cache();

    // Note: my notation vaguely follows that of
    // [wikipedia](http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method).
    // I use the Polak-Ribiere method, with automatic direction reset.
    // Note that we could save some memory by using Fletcher-Reeves, and
    // it seems worth implementing that as an option for
    // memory-constrained problems (then we wouldn't need to store oldgrad).
    if (v >= min_details) {
      printf("\t\tnorm g is %g\n", g.norm());
      printf("\t\toldgrad size is %d\n", oldgrad.get_size());
      printf("\t\tnorm(g - oldgrad) is %g\n", (g-oldgrad).norm());
    }
    double beta = g.dot(g - oldgrad)/oldgradsqr;
    oldgrad = g;
    if (beta < 0 || beta != beta || oldgradsqr == 0) beta = 0;
    oldgradsqr = oldgrad.dot(oldgrad);
    direction = g + beta*direction;
    gdotd = oldgrad.dot(direction);
    if (gdotd < 0) {
      direction = 1.0*oldgrad; // If our direction is uphill, reset to gradient.
      if (v >= verbose) printf("reset to gradient...\n");
      gdotd = oldgrad.dot(direction);
    }
  }

  {
    // Now we will do the line minimization... this is a bit
    // complicated.  We want to use as few steps as possible, but also
    // want to make sure we improve the energy at least a little bit.
  
    const double slope = gdotd;
    const double E0 = energy(quietest(v));
    if (v >= min_details) {
      printf("\t\tQuad: E0 = %25.15g", E0);
      fflush(stdout);
    }
    if (isnan(E0)) {
      // There is no point continuing, since we're starting with a NaN!
      // So we may as well quit here.
      if (v >= verbose) {
        printf(" which is a NaN, so I'm quitting early.\n");
        fflush(stdout);
      }
      return false;
    }
    if (isinf(E0)) {
      // There is no point continuing, since we've got an infinite result.  :(
      // So we may as well quit here.
      if (v >= verbose) {
        printf(" which is infinite, so I'm quitting early.\n");
        fflush(stdout);
      }
      return false;
    }
    if (isnan(slope)) {
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
        printf(" which means we've arrived at the minimum!\n");
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
    *x += step1*direction; // Step away a bit...
    double Etried = -137;
    do {
      step1 *= 0.5;
      *x -= step1*direction; // and then step back a bit less...
      invalidate_cache();
      if (energy(quietest(v)) == Etried) {
        if (v >= verbose) {
          printf("\tThis is silly in QuadraticLineMinimizerType::improve_energy: %g (%g vs %g)\n",
                 step1, energy(quietest(v)), Etried);
          //Grid foo(gd, *x);
          //f.run_finite_difference_test("In QuadraticLineMinimizerType", kT, foo, &direction);
          //fflush(stdout);
        }
        break;
      }
      Etried = energy(quietest(v));
    } while (energy(quietest(v)) > E0 || isnan(energy(quietest(v))));
  
    const double E1 = energy(quietest(v));

    // Do a parabolic extrapolation with E0, E1, and gd to get curvature
    // and new step
    const double curvature = 2.0*(E1-E0-step1*slope)/(step1*step1);
    double step2 = -slope/curvature;
    if (v >= min_details) {
      printf("   E1 = %25.15g\n", E1);
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
        *x -= step1*direction;
        if (v >= verbose) printf("\t\tQuadratic linmin not working properly!!!\n");
        //assert(!"Quadratic linmin not working well!!!\n");
        return false;
      } else {
        // It looks like we aren't moving far enough, so let's try
        // progressively doubling how far we're going.
        double Ebest = E0;
        while (energy(quietest(v)) <= Ebest) {
          Ebest = energy(quietest(v));
          *x += step1*direction;
          invalidate_cache();
          step1 *= 2;
        }
        step1 *= 0.5;
        *x -= step1*direction;
        step = step1;
      }
    } else if (E1 == E0) {
      if (v >= verbose) {
        printf("\t\tNo change in energy... step1 is %g, direction has mag %g pos is %g\n",
               step1, direction.norm(), x->norm());
        fflush(stdout);
      }
      do {
        invalidate_cache();
        *x -= step1*direction;
        step1 *= 2;
        *x += step1*direction;
      } while (energy(quietest(v)) == E0);
      if (energy(quietest(v)) > E0) {
        *x -= step1*direction;
        invalidate_cache();
        step1 = 0;
        printf("\t\tQuad: failed to find any improvement!  :(\n");
      } else if (isnan(energy(quietest(v)))) {
        printf("\t\tQuad: found NaN with smallest possible step!!!\n");
        *x -= step1*direction;
        invalidate_cache();
        step1 = 0;
      }
    } else {
      step = step2; // output the stepsize for later reuse.
      invalidate_cache();
      *x += (step2-step1)*direction; // and move to the expected minimum.
      if (v >= verbose) {
        printf("\t\tQuad: step1 = %14.7g  step2 = %14.7g\n", step1, step2);
        printf("\t\tQuad: E2 = %25.15g\n", energy(quietest(v)));
        fflush(stdout);
      }
      
      // Check that the energy did indeed drop!  FIXME: this may do
      // extra energy calculations, since it's not necessarily shared
      // with the driver routine!
      if (better(E1, energy(quietest(v))) && better(E1, E0)) {
        // The first try was better, so let's go with that one!
        if (v >= verbose) printf("\t\tGoing back to the first try...\n");
        invalidate_cache();
        *x -= (step2-step1)*direction;
      } else if (energy(quietest(v)) > E0) {
        const double E2 = energy(quietest(v));
        invalidate_cache();
        *x -= step2*direction;
        if (v >= verbose) {
          printf("\t\tQuadratic linmin not working well!!! (E2 = %14.7g, E0 = %14.7g)\n",
                 E2, energy(quietest(v)));
          fflush(stdout);
        }
        return false;
        //assert(!"Quadratic linmin not working well!!!\n");
      }
    }
  }
  // Finished with the line minimization!

  if (v >= verbose) {
    //lm->print_info();
    if (v >= min_details) printf("\t\tgrad*direction = %g\n", grad().dot(direction)/gdotd);
    print_info("");
  }
  return (energy(quietest(v)) < E0);
}

void Minimizer::print_info(const char *prefix) const {
  printf("%s==== Iteration %d ====\n", prefix, iter);
  printf("%sEnergy: %g\n\n", prefix, energy(silent));
}
