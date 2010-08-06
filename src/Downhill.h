// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"

class Downhill : public Minimizer {
private:
  double nu;
public:
  Downhill(Functional f, VectorXd *data, double viscosity=0.1)
    : Minimizer(f, data), nu(viscosity) {}

  bool improve_energy(bool verbose = false);
  void print_info(int iter, const char *prefix="") const;
};

class PreconditionedDownhill : public Minimizer {
private:
  double nu;
public:
  PreconditionedDownhill(Functional f, VectorXd *data, double viscosity=0.1)
    : Minimizer(f, data), nu(viscosity) {}

  bool improve_energy(bool verbose = false);
  void print_info(int iter, const char *prefix="") const;
};
