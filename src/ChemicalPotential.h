// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"

class ChemicalPotential : public Functional {
public:
  ChemicalPotential(const GridDescription &g, double chemical_potential)
    : gd(g), mu(chemical_potential) {}

  // To implement a functional, you need to provide both an energy
  // method and a gradient method.
  double operator()(const VectorXd &data) const;
  // If the second pointer is nonzero, you need to also output a
  // preconditioned gradient.
  void grad(const VectorXd &data,
            VectorXd *, VectorXd *pgrad = 0) const;
private:
  GridDescription gd;
  double mu; // the chemical potential
};
