// -*- mode: C++; -*-

#pragma once

#include "lattice.h"

class GridDescription {
public:
  explicit GridDescription(Lattice lat, int nx, int ny, int nz);
  GridDescription(const GridDescription &x);

  Lattice Lat, fineLat;
  int Nx, Ny, Nz, NyNz, NxNyNz;
  double dx, dy, dz;
};
