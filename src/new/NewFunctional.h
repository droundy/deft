// -*- mode: C++; -*-

#pragma once

#include "Vector.h"
#include "Verbosity.h"

struct EnergyGradAndPrecond {
  double energy;
  Vector grad, precond;
};

// A NewFunctional is a very simple type: it defines an energy function
// and a gradient function.  Each of these takes a Vector as input,
// and each takes some additional input to determine its behavior.

// The gradient is slightly more complicated, because it takes an
// optional input that determines which degrees of freedom to take the
// gradient with respect to.  This allows us to either perform a
// constrained minimization, or to perhaps minimize with regard to
// certain degrees of freedom in sequence (if that were more
// efficient).

// The class "NewFunctional" is an abstract class.  To use a NewFunctional,
// you need to create a new class with NewFunctional as its parent.

class NewFunctional {
public:
  NewFunctional() {
    cached_energy = 0;
    data = Vector();
  }
  NewFunctional(const NewFunctional &o) : cached_energy(o.cached_energy), data(o.data) {}
  void operator=(const NewFunctional &o) {
    cached_energy = o.cached_energy;
    data = o.data;
  }
  void operator+=(const Vector &o) {
    invalidate_cache();
    data += o;
  }
  void operator-=(const Vector &o) {
    invalidate_cache();
    data -= o;
  }
  double operator*(const Vector &o) {
    return data.dot(o);
  }
  Vector operator*(double d) {
    return data*d;
  }

  // Use invalidate_cache any time you manually modify the input to
  // the functional.
  void invalidate_cache() const {
    cached_energy = 0;
  }

  // The following is how you should actually find the energy, since
  // it allows us to avoid recomputing the energy when we haven't
  // changed our input.
  double energy() const {
    if (!cached_energy) cached_energy = true_energy();
    return cached_energy;
  }

  // The following is the "true_energy" function, which must be defined by
  // any real functionals.
  virtual double true_energy() const = 0;

  // The following handles the "homeogeneous" case.
  //virtual double energy_per_volume() const = 0;

  // grad is the gradient.  The input "grad_these allows cases where
  // we may want to compute a gradient with regard to just *some* of
  // the input.  This will most likely happen because we are
  // minimizing under some constraint.
  virtual Vector grad() const = 0;

  // The following handles the "homeogeneous" case.
  //virtual double denergy_per_volume_dx() const = 0;

  // We have one additional method, which computes both the grad and
  // the preconditioned gradient.  This one is optional, and only need
  // be defined if there is some actual preconditioning to do;
  virtual EnergyGradAndPrecond energy_grad_and_precond() const {
    assert(0);
  }
  virtual void printme(const char *) const = 0;
  virtual bool have_preconditioner() const { return false; }

  // The following is a utility method to ensure two different
  // functionals have the same input.  Of course, it won't work if
  // they don't expect the same sort of input, so it's mostly only
  // useful for testing.
  void copy_input_data_for_test(const NewFunctional &o) {
    data = o.data;
  }

  // The following is for testing
  int run_finite_difference_test(const char *testname,
                                 const Vector *direction = 0) const;
protected:
  mutable double cached_energy; // store the energy here rather than
                                // recompute it.  zero means invalid.
  mutable Vector data; // yuck, this shouldn't be mutable, but it's a
                       // workaround for the haskell code which is not
                       // yet set up to create non-cost methods.
};
