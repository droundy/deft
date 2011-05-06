// -*- mode: C++; -*-

#include "LineMinimizer.h"
#include <stdio.h>

class QuadraticLineMinimizerType : public MinimizerInterface {
private:
  const VectorXd direction;
  double *step, slope, alternatestep;
public:
  QuadraticLineMinimizerType(Functional f, const GridDescription &gdin, double kT,
                             VectorXd *data, const VectorXd &dir,
                             double gradDotDirection, double *instep)
    : MinimizerInterface(f, gdin, kT, data), direction(dir), step(instep), slope(gradDotDirection) {
    if (step == 0) {
      // If the caller doesn't specify the stepsize, we'll just use an
      // internal one.
      alternatestep = 0.1;
      step = &alternatestep;
    }
  }

  bool improve_energy(bool verbose = false);
  void print_info(int iter) const;
};

bool QuadraticLineMinimizerType::improve_energy(bool verbose) {
  //if (verbose) printf("\t\tI am running QuadraticLineMinimizerType::improve_energy with verbose==%d\n", verbose);
  //fflush(stdout);
  // FIXME: The following probably double-computes the energy!
  const double E0 = energy();
  if (verbose) {
    printf("\t\tQuad: E0 = %25.15g", E0);
    fflush(stdout);
  }
  if (isnan(E0)) {
    // There is no point continuing, since we're starting with a NaN!
    // So we may as well quit here.
    if (verbose) {
      printf(" which is a NaN, so I'm quitting early.\n");
      fflush(stdout);
    }
    return false;
  }
  if (isnan(slope)) {
    // The slope here is a NaN, so there is no point continuing!
    // So we may as well quit here.
    if (verbose) {
      printf(", but the slope is a NaN, so I'm quitting early.\n");
      fflush(stdout);
    }
    return false;
  }

  if (slope*(*step) > 0) {
    if (verbose) printf("Swapping sign of step with slope %g...\n", slope);
    *step *= -1;
  }
  // Here we're going to keep halving step1 until it's small enough
  // that the energy drops a tad.  This hopefully will give us a
  // reasonable starting guess.
  double step1 = *step;
  *x += step1*direction; // Step away a bit...
  double Etried = -137;
  do {
    step1 *= 0.5;
    *x -= step1*direction; // and then step back a bit less...
    invalidate_cache();
    if (energy() == Etried) {
      if (verbose) {
        printf("\tThis is silly in QuadraticLineMinimizerType::improve_energy: %g (%g vs %g)\n",
               step1, energy(), Etried);
        Grid foo(gd, *x);
        f.run_finite_difference_test("In QuadraticLineMinimizerType", kT, foo, &direction);
        fflush(stdout);
      }
      break;
    }
    Etried = energy();
  } while (energy() > E0 || isnan(energy()));
  
  const double E1 = energy();

  // Do a parabolic extrapolation with E0, E1, and gd to get curvature
  // and new step
  const double curvature = 2.0*(E1-E0-step1*slope)/(step1*step1);
  double step2 = -slope/curvature;
  if (verbose) {
    printf("   E1 = %25.15g\n", E1);
    printf("\t\tQuad: slope = %14.7g  curvature = %14.7g\n", slope, curvature);
    fflush(stdout);
  }
 
  if (curvature <= 0.0) {
    if (verbose) {
      printf("\t\tCurvature has wrong sign... %g\n", curvature);
      fflush(stdout);
    }
    if (better(E0, E1)) {
      // It doesn't look to be working, so let's panic!
      invalidate_cache();
      *x -= step1*direction;
      if (verbose) printf("\t\tQuadratic linmin not working properly!!!\n");
      //assert(!"Quadratic linmin not working well!!!\n");
      return false;
    } else {
      // It looks like we aren't moving far enough, so let's try
      // progressively doubling how far we're going.
      double Ebest = E0;
      while (energy() <= Ebest) {
        Ebest = energy();
        *x += step1*direction;
        invalidate_cache();
        step1 *= 2;
      }
      step1 *= 0.5;
      *x -= step1*direction;
      *step = step1;
    }
  } else if (E1 == E0) {
    if (verbose) {
      printf("\t\tNo change in energy...\n");
      fflush(stdout);
    }
    do {
      invalidate_cache();
      *x += step1*direction;
      step1 *= 2;
      //printf("From step1 of %10g I get delta E of %g\n", step1, energy()-E0);
    } while (energy() == E0);
    if (energy() > E0) {
      *x -= step1*direction;
      invalidate_cache();
      step1 = 0;
      printf("\t\tQuad: failed to find any improvement!  :(\n");
    }
  } else {
    *step = step2; // output the stepsize for later reuse.
    invalidate_cache();
    *x += (step2-step1)*direction; // and move to the expected minimum.
    if (verbose) {
      printf("\t\tQuad: step1 = %14.7g  step2 = %14.7g\n", step1, step2);
      printf("\t\tQuad: E2 = %25.15g\n", energy());
      fflush(stdout);
    }

    // Check that the energy did indeed drop!  FIXME: this may do
    // extra energy calculations, since it's not necessarily shared
    // with the driver routine!
    if (better(E1, energy()) && better(E1, E0)) {
      // The first try was better, so let's go with that one!
      if (verbose) printf("\t\tGoing back to the first try...\n");
      invalidate_cache();
      *x -= (step2-step1)*direction;
    } else if (energy() > E0) {
      const double E2 = energy();
      invalidate_cache();
      *x -= step2*direction;
      if (verbose) {
        printf("\t\tQuadratic linmin not working well!!! (E2 = %14.7g, E0 = %14.7g)\n", E2, energy());
        fflush(stdout);
      }
      return false;
      //assert(!"Quadratic linmin not working well!!!\n");
    }
  }
  // We always return false, because it's never recommended to call
  // improve_energy again with this algorithm.
  return false;
}

void QuadraticLineMinimizerType::print_info(int iter) const {
  printf("\t\t==============================\n");
  printf("\t\tLine minimization iteration %2d\n", iter);
  printf("\t\t==============================\n");
}



Minimizer QuadraticLineMinimizer(Functional f, const GridDescription &gd, double kT, VectorXd *data,
                                  const VectorXd &direction, double gradDotDirection, double *step) {
  //if (gradDotDirection > 0) {
  //  printf("The slope is backwards!!!\n");
  //  assert(gradDotDirection < 0);
  //}
  return Minimizer(new QuadraticLineMinimizerType(f, gd, kT, data, direction, gradDotDirection, step));
}
