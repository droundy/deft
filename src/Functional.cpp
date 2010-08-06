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

void FunctionalInterface::print_summary(const char *, const VectorXd &) const {
  // Don't print anything at all by default!
}

void FunctionalInterface::print_iteration(const char *prefix, const VectorXd &d, int iter) const {
  printf("%s==============\n", prefix);
  printf("%sIteration %4d\n", prefix, iter);
  printf("%s==============\n", prefix);
  print_summary(prefix, d);
}

bool FunctionalInterface::run_finite_difference_test(const char *testname,
                                                     const VectorXd &x,
                                                     const VectorXd *direction) const {
  printf("\nRunning finite difference test on %s:\n", testname);

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

  if (min < 1e-3 && best_ratio_error < 1e-7) {
    printf("Passed on basis of reasonable scaling (%g) and accuracy (%g).\n",
           min, best_ratio_error);
    return false;
  }
  if (best_ratio_error < 1e-10) {
    printf("Passed on basis of a gradient accuracy of (%g) (with scaling %g).\n",
           best_ratio_error, min);
    return false;
  }
  if (min < 1e-5 && best_ratio_error < 0.01) {
    printf("Passed on basis of seriously nice scaling %g (and low accuracy %g).\n",
           min, best_ratio_error);
    return false;
  }
  printf("Failed with scaling ratio of %g and gradient accuracy of %g\n",
         min, best_ratio_error);
  return true;
}

class FunctionalComposition : public FunctionalInterface {
public:
  FunctionalComposition(const Functional &f, const FieldFunctional &g)
    : f1(f), f2(g) {};

  double energy(const VectorXd &data) const {
    return f1(f2(data));
  }
  void grad(const VectorXd &data, VectorXd *g, VectorXd *pgrad = 0) const {
    VectorXd d1(f2(data)), g1(d1);
    g1.setZero();
    if (pgrad) {
      VectorXd pg1(g1);
      // First compute the gradient of f1(d1)...
      f1.grad(d1, &g1, &pg1);
      // ... then use the chain rule to compute the gradient of f2(f1(data)).
      // FIXME:  I should be able to use pg1 here!!!
      f2.grad(data, g1, g, pgrad);
    } else {
      // First compute the gradient of f1(d1)...
      f1.grad(d1, &g1);
      // ... then use the chain rule to compute the gradient of f2(f1(data)).
      f2.grad(data, g1, g, pgrad);
    }
  }
  void print_summary(const char *prefix, const VectorXd &data) const {
    f1.print_summary(prefix, f2(data));
  }
private:
  const Functional f1;
  const FieldFunctional f2;
};

Functional compose(const Functional &a, const FieldFunctional &b) {
  return Functional(new FunctionalComposition(a, b));
}

// The sum of two functionals...

class FunctionalSum : public FunctionalInterface {
public:
  FunctionalSum(const Functional &f, const Functional &g) : f1(f), f2(g) {};

  double energy(const VectorXd &data) const {
    return f1(data) + f2(data);
  }
  void grad(const VectorXd &data, VectorXd *g, VectorXd *pgrad) const {
    f1.grad(data, g, pgrad);
    f2.grad(data, g, pgrad);
  }
  void print_summary(const char *prefix, const VectorXd &data) const {
    f1.print_summary(prefix, data);
    f2.print_summary(prefix, data);
  }
private:
  const Functional f1, f2;
};

Functional operator+(const Functional &a, const Functional &b) {
  return Functional(new FunctionalSum(a, b));
}
