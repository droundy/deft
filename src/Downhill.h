// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"

class Downhill : public Minimizer {
private:
  double nu;
public:
  Downhill(counted_ptr<Functional> f, VectorXd *data, double viscosity=0.1)
    : Minimizer(f, data), nu(viscosity) {}

  bool improve_energy(bool verbose = false);
  void print_info(int iter) const;
};
