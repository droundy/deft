// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"

class ExternalPotential : public Functional {
public:
  ExternalPotential(const Grid &V)
    : Vdvolume(V*V.description().Lat.volume()/V.description().NxNyNz) {}

  // To implement a functional, you need to provide both an energy
  // method and a gradient method.
  double operator()(const VectorXd &data) const;
  // If the second pointer is nonzero, you need to also output a
  // preconditioned gradient.
  void grad(const VectorXd &data,
            VectorXd *, VectorXd *pgrad = 0) const;

  // You may optionally define a print_summary method, which would
  // print something interesting to the screen.
  void print_summary() const;
private:
  VectorXd Vdvolume; // the potential times the differential volume element dV.
};
