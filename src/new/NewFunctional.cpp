#include "NewFunctional.h"
#include <math.h>
#include <float.h>

int NewFunctional::run_finite_difference_test(const char *testname, const Vector *direction,
                                              double stepsize) const {
  printf("\nRunning finite difference test on %s:\n", testname);

  int retval = 0;

  Vector my_grad = grad();

  double Eold = energy();
  Vector my_direction(my_grad.get_size());
  if (direction) my_direction = *direction;
  else my_direction = my_grad;
  my_direction /= my_direction.norm();

  const double lderiv = my_direction.dot(my_grad);
  if (lderiv == 0.0) {
    printf("FAIL: Gradient is zero...\n");
    retval++;
  }

  // We try to choose as small values for epsilon as are consistent with
  // getting meaningful results:
  const double sigfigs_best = 1e-15;
  const double epsilon_max = 1e-15*fabs(Eold)/fabs(lderiv)/sigfigs_best;
  int min_p = (int) -floor(log10(epsilon_max) + 0.5);
  if (stepsize) {
    const int old_min_p = min_p;
    min_p = -floor(log10(stepsize*my_grad.norm()) + 0.5);
    printf("user-specified delta of 1e%d instead of 1e%d from directional grad %g and energy %g\n",
           -min_p, -old_min_p, lderiv, Eold);
  } else {
    printf("choosing delta of 1e%d based on directional grad %g and energy %g\n",
           -min_p, lderiv, Eold);
  }
  const int max_p = min_p + 7;
  double min = 1e300, best_ratio_error = 1.0;
  Vector grads(max_p + 1 - min_p);
  Vector diff_grads(max_p - min_p);
  // Take steps in the grad direction.
  for(int p=min_p; p <= max_p; p++) {
    const double eps_ratio = 10.0;
    const double epsilon = pow(eps_ratio, -p);
    // The following is a little wasteful of memory...
    data += epsilon*my_direction;
    const double Eplus=energy();
    data -= 2*epsilon*my_direction;
    const double Eminus=energy();
    data += epsilon*my_direction;

    printf("Eplus %g  Eminus %g\n", Eplus, Eminus);
    grads[p-min_p] = (Eplus-Eminus)/(2*epsilon);
    printf("    eps^2 = %25.16f deltaE %.12g\n",
           epsilon*epsilon*pow(eps_ratio, 2.0*min_p), Eplus-Eminus);
    printf("FD   Ratio: %25.16f (grad ~ %g)\n", grads[p-min_p]/lderiv,
           grads[p-min_p]);
    printf("FD sigfigs: %25.16f\n", 1e-15*fabs(Eold/(Eplus-Eminus)));
    fflush(stdout);

    if (fabs(grads[p-min_p]/lderiv - 1.0) < best_ratio_error)
      best_ratio_error = fabs(grads[p-min_p]/lderiv - 1.0);
    if (p > min_p) {
      diff_grads[p-min_p-1]=(grads[p-min_p]-grads[p-min_p-1]);
      double diff = (grads[p-min_p-1]-lderiv)*
        (-1+1./(eps_ratio*eps_ratio))/diff_grads[p-min_p-1];
      if ( min > fabs(diff-1) ) min = fabs(diff-1);
      //dft_log("FD diff: %25.16lf\n\n", diff);
    }
  }

  if (min < 1e-3 && best_ratio_error < 1e-7) {
    printf("Passed on basis of reasonable scaling (%g) and accuracy (%g).\n",
           min, best_ratio_error);
  } else if (best_ratio_error < 1e-10) {
    printf("Passed on basis of a gradient accuracy of (%g) (with scaling %g).\n",
           best_ratio_error, min);
  } else if (min < 1e-5 && best_ratio_error < 0.01) {
    printf("Passed on basis of seriously nice scaling %g (and low accuracy %g).\n",
           min, best_ratio_error);
  } else {
    printf("FAIL: Failed with scaling ratio of %g and gradient accuracy of %g\n",
           min, best_ratio_error);
    retval++;
  }

  return retval;
}
