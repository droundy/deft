// -*- mode: C++; -*-

#pragma once

#include "Vector.h"
#include "Bools.h"

struct EnergyGradAndPrecond {
  double energy;
  Vector grad, precond;
};

// A Functional is a very simple type: it defines an energy function
// and a gradient function.  Each of these takes a Vector as input,
// and each takes some additional input to determine its behavior.

// The gradient is slightly more complicated, because it takes an
// optional input that determines which degrees of freedom to take the
// gradient with respect to.  This allows us to either perform a
// constrained minimization, or to perhaps minimize with regard to
// certain degrees of freedom in sequence (if that were more
// efficient).

// The class "Functional" is an abstract class.  To use a Functional,
// you need to create a new class with Functional as its parent.

class Functional {
public:
  // The following is the "energy" function, which must be defined by
  // any real functionals.
  virtual double energy(const Vector &x, bool verbose = false) const = 0;

  // grad is the gradient.  The input "grad_these allows cases where
  // we may want to compute a gradient with regard to just *some* of
  // the input.  This will most likely happen because we are
  // minimizing under some constraint.
  virtual Vector grad(const Vector &x, const Bools *grad_these = 0) const = 0;

  // We have one additional method, which computes both the grad and
  // the preconditioned gradient.  This one is optional, and only need
  // be defined if there is some actual preconditioning to do;
  virtual EnergyGradAndPrecond energy_grad_and_precond(const Vector &x,
                                                       bool verbose,
                                                       const Bools *grad_these) const {
    assert(0);
  }
  virtual bool have_preconditioner() const { return false; }
};
