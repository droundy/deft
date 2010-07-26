// -*- mode: C++; -*-

#pragma once

#include <eigen2/Eigen/Core>

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class Functional {
public:
  Functional(VectorXd *data) : mydata(data) {}

  VectorXd &data() { return *mydata; }
  const VectorXd &data() const { return *mydata; }

  // To implement a functional, you need to provide both an energy
  // method and a gradient method.
  virtual double energy() const = 0;
  // If the second pointer is nonzero, you need to also output a
  // preconditioned gradient.
  virtual void grad(VectorXd *, VectorXd *pgrad = 0) const = 0;

  // You may optionally define a print_iteration_summary method, which
  // would print something interesting to the screen.
  virtual void  print_summary() const;

  // The following utility methods do not need to be overridden.
  void print_iteration(int iter) const;
  bool run_finite_difference_test(const char *testname,
                                  const VectorXd *direction = 0);
private:
  VectorXd *mydata;
};
