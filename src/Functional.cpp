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

#include "Functional.h"
#include "handymath.h"
#include <stdio.h>
#include <math.h>

void Functional::print_summary(const VectorXd &) const {
  printf("hello world\n");
}

void Functional::print_iteration(const VectorXd &d, int iter) const {
  printf("Iteration number %d:\n", iter);
  print_summary(d);
}

bool Functional::run_finite_difference_test(const char *testname,
                                            const VectorXd &x,
                                            const VectorXd *direction) {
  printf("Running finite difference test on %s:\n", testname);

  VectorXd my_grad(x);
  double Eold = (*this)(x);
  my_grad.setZero();
  grad(x, &my_grad);
  VectorXd my_direction(my_grad);
  if (direction) my_direction = *direction;

  const double lderiv = my_direction.dot(my_grad);
  if (lderiv == 0.0) {
    printf("Gradient is zero...\n");
    return true;
  }
  {
    // compare the gradient computed by the two functions grad() and
    // grad_and_pgrad()
    VectorXd my_pgrad(x);
    my_grad.setZero();
    my_pgrad.setZero();
    grad(x, &my_grad, &my_pgrad);
    const double lderiv_new = my_direction.dot(my_grad);
    if (fabs(lderiv_new/lderiv - 1) > 1e-12) {
      printf("\n*** WARNING!!! INCONSISTENT GRADIENTS! ***\n");
      printf("Different gradient in energy_and_grad() and grad_and_pgrad()\n");
      printf("Fractional error is %g\n\n", lderiv_new/lderiv - 1);
      return false;
    }
  }

  // We try to choose as small values for epsilon as are consistent with
  // getting meaningful results:
  const double sigfigs_best = 1e-15;
  const double epsilon_max = 1e-15*fabs(Eold)/fabs(lderiv)/sigfigs_best;
  const int min_p = (int) -floor(log10(epsilon_max) + 0.5);
  const int max_p = min_p + 7;
  printf("choosing delta of grad*1e%d based on grad %g and energy %g\n",
         -min_p, lderiv, Eold);
  double min = 1e300, best_ratio_error = 1.0;
  VectorXd grads(max_p + 1 - min_p);
  VectorXd diff_grads(max_p - min_p);
  // Take steps in the grad direction.
  for(int p=min_p; p <= max_p; p++) {
    const double eps_ratio = 10.0;
    const double epsilon = pow(eps_ratio, -p);
    // The following is a little wasteful of memory...
    const double Eplus= (*this)(x + epsilon*my_direction);
    const double Eminus=(*this)(x - epsilon*my_direction);

    grads[p-min_p] = (Eplus-Eminus)/(2*epsilon);
    printf("    eps^2 = %25.16f deltaE %.12g\n",
           epsilon*epsilon*pow(eps_ratio, 2.0*min_p), Eplus-Eminus);
    printf("FD   Ratio: %25.16f (grad ~ %g)\n", grads[p-min_p]/lderiv,
           grads[p-min_p]);
    printf("FD sigfigs: %25.16f\n", 1e-15*fabs(Eold/(Eplus-Eminus)));

    if (fabs(grads[p-min_p]/lderiv - 1.0) < best_ratio_error)
      best_ratio_error = fabs(grads[p-min_p]/lderiv - 1.0);
    if (p > min_p) {
      diff_grads[p-min_p-1]=(grads[p-min_p]-grads[p-min_p-1]);
      double diff = (grads[p-min_p-1]-lderiv)*
        (-1+1./sqr(eps_ratio))/diff_grads[p-min_p-1];
      if ( min > fabs(diff-1) ) min = fabs(diff-1);
      //dft_log("FD diff: %25.16lf\n\n", diff);
    }
  }

  //dft_log("The min is %g and best_ratio_error is %g\n",
  //        min, best_ratio_error);
  return (min < 1e-3) || (best_ratio_error < 1e-9);

}
