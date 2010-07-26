// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"

class IdealGas : public Functional {
public:
  IdealGas(Grid *data, double temp)
    : Functional(data), n(*data), T(temp) {}

  // To implement a functional, you need to provide both an energy
  // method and a gradient method.
  double energy() const;
  // If the second pointer is nonzero, you need to also output a
  // preconditioned gradient.
  void grad(VectorXd *, VectorXd *pgrad = 0) const;

  // You may optionally define a print_iteration_summary method, which
  // would print something interesting to the screen.
  void  print_summary() const;
private:
  Grid &n;
  double T; // temperature
};
