// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"

class Downhill : public Minimizer {
protected:
  double nu, orig_nu;
public:
  Downhill(Functional f, const GridDescription &gdin, VectorXd *data, double viscosity=0.1)
    : Minimizer(f, gdin, data), nu(viscosity), orig_nu(viscosity) {}
  void minimize(Functional newf, const GridDescription &gdnew, VectorXd *newx = 0) {
    nu = orig_nu;
    Minimizer::minimize(newf, gdnew, newx);
  }

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};

class PreconditionedDownhill : public Downhill {
public:
  PreconditionedDownhill(Functional f, const GridDescription &gdin, VectorXd *data, double viscosity=0.1)
    : Downhill(f, gdin, data, viscosity) {}

  bool improve_energy(bool verbose = false);
};
