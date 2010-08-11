// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"

class Downhill : public Minimizer {
private:
  double nu, orig_nu;
public:
  Downhill(Functional f, VectorXd *data, double viscosity=0.1)
    : Minimizer(f, data), nu(viscosity), orig_nu(viscosity) {}
  void minimize(Functional newf, VectorXd *newx = 0) {
    nu = orig_nu;
    Minimizer::minimize(newf, newx);
  }

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};

class PreconditionedDownhill : public Minimizer {
private:
  double nu, orig_nu;
public:
  PreconditionedDownhill(Functional f, VectorXd *data, double viscosity=0.1)
    : Minimizer(f, data), nu(viscosity), orig_nu(viscosity) {}
  void minimize(Functional newf, VectorXd *newx = 0) {
    nu = orig_nu;
    Minimizer::minimize(newf, newx);
  }

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};
