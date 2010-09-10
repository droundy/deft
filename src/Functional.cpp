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
#include <math.h>

void FunctionalInterface::print_summary(const char *, double) const {
  // Don't print anything at all by default!
}

void Functional::print_summary(const char *prefix) const {
  itsCounter->ptr->print_summary(prefix, itsCounter->last_energy);
  if (itsCounter->name) {
    printf("%s%25s =", prefix, itsCounter->name);
    print_double("", itsCounter->last_energy);
    printf("\n");
  }
}

bool Functional::run_finite_difference_test(const char *testname, const Grid &x,
                                            const VectorXd *direction) const {
  printf("\nRunning finite difference test on %s:\n", testname);
  const GridDescription gd(x.description());

  VectorXd my_grad(x);
  double Eold = (*this)(x);
  my_grad.setZero();
  grad(x, &my_grad);
  VectorXd my_direction(my_grad);
  if (direction) my_direction = *direction;
  my_direction /= my_direction.norm();

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
  printf("choosing delta of 1e%d based on directional grad %g and energy %g\n",
         -min_p, lderiv, Eold);
  double min = 1e300, best_ratio_error = 1.0;
  VectorXd grads(max_p + 1 - min_p);
  VectorXd diff_grads(max_p - min_p);
  // Take steps in the grad direction.
  for(int p=min_p; p <= max_p; p++) {
    const double eps_ratio = 10.0;
    const double epsilon = pow(eps_ratio, -p);
    // The following is a little wasteful of memory...
    const double Eplus= (*this)(gd, x + epsilon*my_direction);
    const double Eminus=(*this)(gd, x - epsilon*my_direction);

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
  printf("FAIL: Failed with scaling ratio of %g and gradient accuracy of %g\n",
         min, best_ratio_error);
  return true;
}

void Functional::print_iteration(const char *prefix, int iter) const {
  printf("%s==============\n", prefix);
  printf("%sIteration %4d\n", prefix, iter);
  printf("%s==============\n", prefix);
  print_summary(prefix);
}

class FunctionalComposition : public FunctionalInterface {
public:
  FunctionalComposition(const Functional &f, const FieldFunctional &g)
    : f1(f), f2(g) {};

  double energy(const GridDescription &gd, const VectorXd &data) const {
    return f1(gd, f2(gd, data));
  }
  double energy(double data) const {
    return f1(f2(data));
  }
  double grad(double data) const {
    return f1.grad(f2(data))*f2.grad(data);
  }
  void grad(const GridDescription &gd, const VectorXd &data, VectorXd *g, VectorXd *pgrad = 0) const {
    VectorXd d1(f2(gd, data)), g1(d1);
    g1.setZero();
    if (pgrad) {
      VectorXd pg1(g1);
      // First compute the gradient of f1(d1)...
      f1.grad(gd, d1, &g1, &pg1);
      // ... then use the chain rule to compute the gradient of f2(f1(data)).
      // FIXME:  I should be able to use pg1 here!!!
      f2.grad(gd, data, g1, g, pgrad);
    } else {
      // First compute the gradient of f1(d1)...
      f1.grad(gd, d1, &g1);
      // ... then use the chain rule to compute the gradient of f2(f1(data)).
      f2.grad(gd, data, g1, g, 0);
    }
  }
  void print_summary(const char *prefix, double) const {
    f1.print_summary(prefix);
  }
private:
  const Functional f1;
  const FieldFunctional f2;
};

Functional Functional::operator()(const FieldFunctional &f) const {
  return Functional(new FunctionalComposition(*this, f), f.get_name());
}

// The sum of two functionals...

class FunctionalSum : public FunctionalInterface {
public:
  FunctionalSum(const Functional &f, const Functional &g) : f1(f), f2(g) {};
  ~FunctionalSum() {}

  double energy(double data) const {
    return f1(data) + f2(data);
  }
  double grad(double data) const {
    return f1.grad(data) + f2.grad(data);
  }
  double energy(const GridDescription &gd, const VectorXd &data) const {
    return f1(gd, data) + f2(gd, data);
  }
  void grad(const GridDescription &gd, const VectorXd &data, VectorXd *g, VectorXd *pgrad) const {
    f1.grad(gd, data, g, pgrad);
    f2.grad(gd, data, g, pgrad);
  }
  void print_summary(const char *prefix, double) const {
    f1.print_summary(prefix);
    f2.print_summary(prefix);
  }
private:
  const Functional f1, f2;
};

Functional operator+(const Functional &a, const Functional &b) {
  return Functional(new FunctionalSum(a, b));
}


// The product of two functionals...

class FunctionalProduct : public FunctionalInterface {
public:
  FunctionalProduct(const Functional &f, const Functional &g) : f1(f), f2(g) {};

  double energy(double data) const {
    return f1(data)*f2(data);
  }
  double grad(double data) const {
    return f1(data)*f2.grad(data) + f1.grad(data)*f2(data);
  }
  double energy(const GridDescription &gd, const VectorXd &data) const {
    return f1(gd, data)*f2(gd, data);
  }
  void grad(const GridDescription &gd, const VectorXd &data, VectorXd *g, VectorXd *pgrad) const {
    Grid tempgrad(gd);
    tempgrad.setZero();
    if (pgrad) {
      Grid tempp(gd);
      tempp.setZero();
      f1.grad(gd, data, &tempgrad, &tempp);
      double f2data = f2(gd, data);
      *g += f2data*tempgrad;
      *pgrad += f2data*tempp;
      tempgrad.setZero();
      tempp.setZero();
      f2.grad(gd, data, &tempgrad, &tempp);
      double f1data = f1(gd, data);
      *g += f1data*tempgrad;
      *pgrad += f1data*tempp;
    } else {
      f1.grad(gd, data, &tempgrad, 0);
      *g += f2(gd, data)*tempgrad;
      tempgrad.setZero();
      f2.grad(gd, data, &tempgrad, 0);
      *g += f1(gd, data)*tempgrad;
    }
  }
  void print_summary(const char *prefix, double) const {
    f1.print_summary(prefix);
    f2.print_summary(prefix);
  }
private:
  const Functional f1, f2;
};

Functional operator*(const Functional &a, const Functional &b) {
  return Functional(new FunctionalProduct(a, b));
}


// The constant times a functional...

class ScalarProduct : public FunctionalInterface {
public:
  ScalarProduct(double x, const Functional &y) : a(x), f(y) {};

  double energy(double data) const {
    return a*f(data);
  }
  double grad(double data) const {
    return a*f.grad(data);
  }
  double energy(const GridDescription &gd, const VectorXd &data) const {
    return a*f(gd, data);
  }
  void grad(const GridDescription &gd, const VectorXd &data, VectorXd *g, VectorXd *pgrad) const {
    Grid tempgrad(gd);
    tempgrad.setZero();
    if (pgrad) {
      Grid tempp(gd);
      tempp.setZero();
      f.grad(gd, data, &tempgrad, &tempp);
      *g += a*tempgrad;
      *pgrad += a*tempp;
    } else {
      f.grad(gd, data, &tempgrad, 0);
      *g += a*tempgrad;
    }
  }
private:
  double a;
  const Functional f;
};

Functional operator*(double a, const Functional &b) {
  return Functional(new ScalarProduct(a, b));
}


// A simple integrating functional...

class Integrate : public FunctionalInterface {
public:
  Integrate(const FieldFunctional &ff) : f(ff) {};

  double energy(double data) const {
    return f(data);
  }
  double grad(double data) const {
    return f.grad(data);
  }
  double energy(const GridDescription &gd, const VectorXd &data) const {
    // This does some extra work to save the energies of each term in
    // the sum.
    VectorXd fdata(f.justMe(gd, data));
    double e = gd.dvolume*fdata.sum();
    f.last_energy = e;
    FieldFunctional *nxt = f.next();
    while (nxt) {
      fdata += nxt->justMe(gd, data);
      double etot = gd.dvolume*fdata.sum();
      nxt->last_energy = etot - e;
      e = etot;
      nxt = nxt->next();
    }
    return e;
  }
  void grad(const GridDescription &gd, const VectorXd &x, VectorXd *g, VectorXd *pgrad) const {
    f.grad(gd, x, gd.dvolume*VectorXd::Ones(gd.NxNyNz), g, pgrad);
  }
  void print_summary(const char *prefix, double) const {
    const FieldFunctional *nxt = &f;
    while (nxt) {
      if (nxt->get_name()) {
        printf("%s%25s =", prefix, nxt->get_name());
        print_double("", nxt->last_energy);
        printf("\n");
      } else {
        printf("%s%25s =", prefix, "UNKNOWN");
        print_double("", nxt->last_energy);
        printf("\n");
      }
      nxt = nxt->next();
    }
  }
private:
  const FieldFunctional f;
};

Functional integrate(const FieldFunctional &f) {
  return Functional(new Integrate(f));
}

// The constant times a functional...

class Constraint : public FunctionalInterface {
public:
  Constraint(const Grid &g, const Functional &y) : constraint(g), f(y) {};

  double energy(double data) const {
    return f(data);
  }
  double grad(double data) const {
    return f.grad(data);
  }
  double energy(const GridDescription &gd, const VectorXd &data) const {
    return f(gd, data);
  }
  void grad(const GridDescription &gd, const VectorXd &data, VectorXd *g, VectorXd *pgrad) const {
    f.grad(gd, data, g, pgrad);
    g->cwise() *= constraint;
    if (pgrad) pgrad->cwise() *= constraint;
  }
  void print_summary(const char *prefix, double) const {
    f.print_summary(prefix);
  }
private:
  const Grid constraint;
  const Functional f;
};

Functional constrain(const Grid &g, Functional f) {
  return Functional(new Constraint(g, f));
}
