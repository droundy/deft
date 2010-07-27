// -*- mode: C++; -*-

#pragma once

#include "FieldFunctional.h"
#include "counted_ptr.h"

class FunctionalComposition;

class Functional {
public:
  // To implement a functional, you need to provide both an energy
  // method and a gradient method.

  // The energy method is operator()
  virtual double operator()(const VectorXd &data) const = 0;

  // The gradient method, however, outputs its output in-place in two
  // VectorXd pointers.  If the second pointer is nonzero, you need to
  // also output a preconditioned gradient.  The two gradients should
  // be accumulated, so that the resulting VectorXd will be the sum of
  // the gradient plus any additional value.  This enables the easy
  // combination of different functionals.
  virtual void grad(const VectorXd &data,
                    VectorXd *, VectorXd *pgrad = 0) const = 0;

  // You may optionally define a print_iteration_summary method, which
  // would print something interesting to the screen.
  virtual void print_summary(const VectorXd &data) const;

  // The following utility methods do not need to be overridden.
  void print_iteration(const VectorXd &data, int iter) const;
  bool run_finite_difference_test(const char *testname,
                                  const VectorXd &data,
                                  const VectorXd *direction = 0);
};

class FunctionalComposition : public Functional {
public:
  FunctionalComposition(const counted_ptr<Functional> f,
                        const counted_ptr<FieldFunctional> g)
    : f1(f), f2(g) {};

  double operator()(const VectorXd &data) const;
  void grad(const VectorXd &data, VectorXd *, VectorXd *pgrad = 0) const;
  void print_summary(const VectorXd &data) const;
private:
  const counted_ptr<Functional> f1;
  const counted_ptr<FieldFunctional> f2;
};

counted_ptr<FunctionalComposition>
    compose(const counted_ptr<Functional> a,
            const counted_ptr<FieldFunctional> b);
