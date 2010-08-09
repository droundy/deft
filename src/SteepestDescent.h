// -*- mode: C++; -*-

#pragma once

#include "LineMinimizer.h"

class SteepestDescent : public Minimizer {
private:
  double step;
  LineMinimizer linmin;
public:
  SteepestDescent(Functional f, VectorXd *data, LineMinimizer lm, double stepsize = 0.1)
    : Minimizer(f, data), step(stepsize), linmin(lm) {}

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};

class PreconditionedSteepestDescent : public Minimizer {
private:
  double step;
  LineMinimizer linmin;
public:
  PreconditionedSteepestDescent(Functional f, VectorXd *data, LineMinimizer lm, double stepsize = 0.1)
    : Minimizer(f, data), step(stepsize), linmin(lm) {}

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};
