// -*- mode: C++; -*-

#pragma once

#include "LineMinimizer.h"

class SteepestDescent : public Minimizer {
private:
  double step, orig_step;
  LineMinimizer linmin;
public:
  SteepestDescent(Functional f, VectorXd *data, LineMinimizer lm, double stepsize = 0.1)
    : Minimizer(f, data), step(stepsize), orig_step(step), linmin(lm) {}
  void minimize(Functional newf, VectorXd *newx = 0) {
    step = orig_step;
    Minimizer::minimize(newf, newx);
  }

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};

class PreconditionedSteepestDescent : public Minimizer {
private:
  double step, orig_step;
  LineMinimizer linmin;
public:
  PreconditionedSteepestDescent(Functional f, VectorXd *data, LineMinimizer lm, double stepsize = 0.1)
    : Minimizer(f, data), step(stepsize), orig_step(step), linmin(lm) {}
  void minimize(Functional newf, VectorXd *newx = 0) {
    step = orig_step;
    Minimizer::minimize(newf, newx);
  }

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};
