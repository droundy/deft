#include "GridDescription.h"
#include "Grid.h"
#include "ReciprocalGrid.h"
#include <fftw3.h>

GridDescription::GridDescription(Lattice lat, int nx, int ny, int nz)
  : Lat(lat), fineLat(Cartesian(lat.a1()/nx), Cartesian(lat.a2()/ny),
                      Cartesian(lat.a3()/nz)) {
  Nx = nx; Ny = ny; Nz = nz;
  NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
  NzOver2 = Nz/2 + 1; NyNzOver2 = Ny*NzOver2; NxNyNzOver2 = Nx*NyNzOver2;
  dx = 1.0/Nx; dy = 1.0/Ny; dz = 1.0/Nz;
  dvolume = fineLat.volume();

  // Make a couple of FFTW plans with FFTW_MEASURE, to speed things up
  // for later...
  Grid rspace(*this);
  ReciprocalGrid kspace(*this);
  fftw_destroy_plan(fftw_plan_dft_r2c_3d(Nx, Ny, Nz, rspace.data(),
                                         (fftw_complex *)kspace.data(),
                                         FFTW_MEASURE));
  fftw_destroy_plan(fftw_plan_dft_c2r_3d(Nx, Ny, Nz,
                                         (fftw_complex *)kspace.data(),
                                         rspace.data(), FFTW_MEASURE));
}

GridDescription::GridDescription(Lattice lat, double delta)
  : Nx(1+int(lat.a1().norm()/delta)),
    Ny(1+int(lat.a2().norm()/delta)),
    Nz(1+int(lat.a3().norm()/delta)),
    Lat(lat), fineLat(Cartesian(lat.a1()/Nx), Cartesian(lat.a2()/Ny),
                      Cartesian(lat.a3()/Nz)) {
  NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
  NzOver2 = Nz/2 + 1; NyNzOver2 = Ny*NzOver2; NxNyNzOver2 = Nx*NyNzOver2;
  dx = 1.0/Nx; dy = 1.0/Ny; dz = 1.0/Nz;
  dvolume = fineLat.volume();

  // Make a couple of FFTW plans with FFTW_MEASURE, to speed things up
  // for later...
  Grid rspace(*this);
  ReciprocalGrid kspace(*this);
  fftw_destroy_plan(fftw_plan_dft_r2c_3d(Nx, Ny, Nz, rspace.data(),
                                         (fftw_complex *)kspace.data(),
                                         FFTW_MEASURE));
  fftw_destroy_plan(fftw_plan_dft_c2r_3d(Nx, Ny, Nz,
                                         (fftw_complex *)kspace.data(),
                                         rspace.data(), FFTW_MEASURE));
}
