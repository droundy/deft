// -*- mode: C++; -*-

#pragma once

#include "LineMinimizer.h"

class SteepestDescent : public Minimizer {
protected:
  double step, orig_step;
  LineMinimizer linmin;
public:
  SteepestDescent(Functional f, const GridDescription &gdin, VectorXd *data, LineMinimizer lm, double stepsize = 0.1)
    : Minimizer(f, gdin, data), step(stepsize), orig_step(step), linmin(lm) {}
  void minimize(Functional newf, const GridDescription &gdnew, VectorXd *newx = 0) {
    step = orig_step;
    Minimizer::minimize(newf, gdnew, newx);
  }

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};

class PreconditionedSteepestDescent : public SteepestDescent {
public:
  PreconditionedSteepestDescent(Functional f, const GridDescription &gdin,
                                VectorXd *data, LineMinimizer lm, double stepsize = 0.1)
    : SteepestDescent(f, gdin, data, lm, stepsize) {}

  bool improve_energy(bool verbose = false);
};
