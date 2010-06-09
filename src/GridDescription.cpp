#include "GridDescription.h"

GridDescription::GridDescription(Lattice lat, int nx, int ny, int nz)
  : Lat(lat), fineLat(Cartesian(lat.a1()/nx), Cartesian(lat.a2()/ny),
                      Cartesian(lat.a3()/nz)) {
  Nx = nx; Ny = ny; Nz = nz;
  NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
  dx = 1.0/Nx; dy = 1.0/Ny; dz = 1.0/Nz;
}

GridDescription::GridDescription(const GridDescription &x)
  : Lat(x.Lat), fineLat(x.fineLat) {
  Nx = x.Nx; Ny = x.Ny; Nz = x.Nz;
  NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
  dx = 1.0/Nx; dy = 1.0/Ny; dz = 1.0/Nz;
}
