// -*- mode: C++; -*-

#pragma once

#include "lattice.h"

class GridDescription {
public:
  explicit GridDescription(Lattice lat, int nx, int ny, int nz);
  // Default copy constructor is just fine!

  Lattice Lat, fineLat;
  int Nx, Ny, Nz, NyNz, NxNyNz, NzOver2, NyNzOver2, NxNyNzOver2;
  double dx, dy, dz;
private:
  void initme();
};
