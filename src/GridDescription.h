// -*- mode: C++; -*-

#pragma once

#include "lattice.h"

#pragma GCC diagnostic push
#if __GCC__ > 5
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include "Eigen/Eigen"
#pragma GCC diagnostic pop


typedef std::complex<double> complex;

class GridDescription {
public:
  explicit GridDescription(Lattice lat, int nx, int ny, int nz);
  explicit GridDescription(Lattice lat, double dx);
  // Default copy constructor is just fine!

  double dx, dy, dz, dvolume;
  int Nx, Ny, Nz, NyNz, NxNyNz, NzOver2, NyNzOver2, NxNyNzOver2;
  Lattice Lat, fineLat;
private:
  void initme();
};

#include "ReciprocalOperators.h"
#include "RealSpaceOperators.h"
