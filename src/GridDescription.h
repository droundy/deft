// -*- mode: C++; -*-

#pragma once

#include "lattice.h"

class GridDescription {
public:
  explicit GridDescription(Lattice lat, int nx, int ny, int nz);
  explicit GridDescription(Lattice lat, double dx);
  // Default copy constructor is just fine!

  double dx, dy, dz;
  int Nx, Ny, Nz, NyNz, NxNyNz, NzOver2, NyNzOver2, NxNyNzOver2;
  Lattice Lat, fineLat;
private:
  void initme();
};
