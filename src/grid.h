#pragma once

#include "lattice.h"

class Grid : public VectorXd {
public:
  explicit Grid(Lattice lat, int nx, int ny, int nz)
    : VectorXd(nx*ny*nz), Lat(lat) {
    Nx = nx; Ny = ny; Nz = nz;
    NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
  }
  Grid(const Grid &x) : VectorXd(x.Nx*x.Ny*x.Nz), Lat(x.Lat) {
    Nx = x.Nx; Ny = x.Ny; Nz = x.Nz;
    NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
  }
  double operator()(int x, int y, int z) const {
    return (*this)[x*NyNz + y*Nz + z];
  }
  double &operator()(int x, int y, int z) {
    return (*this)[x*NyNz + y*Nz + z];
  }
  double operator()(const Cartesian &r) const {
    return (*this)(Lat.toRelative(r));
  }
  double operator()(const Relative &r) const {
    double rx = r(0)*Nx, ry = r(1)*Ny, rz = r(2)*Nz;
    int ix = int(rx), iy = int(ry), iz = int(rz);
    int ixp1 = (ix+1)%Nx, iyp1 = (iy+1)%Ny, izp1 = (iz+1)%Nz;
    double wx = rx-ix, wy = ry-iy, wz = rz-iz;
    return (1-wx)*(1-wy)*(1-wz)*(*this)(ix,iy,iz)
      + wx*(1-wy)*(1-wz)*(*this)(ixp1,iy,iz)
      + (1-wx)*wy*(1-wz)*(*this)(ix,iyp1,iz)
      + (1-wx)*(1-wy)*wz*(*this)(ix,iy,izp1)
      + wx*(1-wy)*wz*(*this)(ixp1,iy,izp1)
      + (1-wx)*wy*wz*(*this)(ix,iyp1,izp1)
      + wx*wy*(1-wz)*(*this)(ixp1,iyp1,iz)
      + wx*wy*wz*(*this)(ixp1,iyp1,izp1);
  }
private:
  Lattice Lat;
  int Nx, Ny, Nz, NyNz, NxNyNz;
};

class ReciprocalGrid : public VectorXd {
public:
  explicit ReciprocalGrid(Lattice lat, int nx, int ny, int nz)
    : VectorXd(nx*ny*nz), Lat(lat) {
    Nx = nx; Ny = ny; Nz = nz;
    NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
  }
  ReciprocalGrid(const ReciprocalGrid &x)
    : VectorXd(x.Nx*x.Ny*x.Nz), Lat(x.Lat) {
    Nx = x.Nx; Ny = x.Ny; Nz = x.Nz;
    NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
  }
  double operator()(int x, int y, int z) const {
    return (*this)[x*NyNz + y*Nz + z];
  }
  double &operator()(int x, int y, int z) {
    return (*this)[x*NyNz + y*Nz + z];
  }
  double operator()(const Reciprocal &r) const {
    return (*this)(Lat.toRelativeReciprocal(r));
  }
  double operator()(const RelativeReciprocal &r) const {
    double rx = r(0)*Nx, ry = r(1)*Ny, rz = r(2)*Nz;
    int ix = int(rx), iy = int(ry), iz = int(rz);
    int ixp1 = (ix+1)%Nx, iyp1 = (iy+1)%Ny, izp1 = (iz+1)%Nz;
    double wx = rx-ix, wy = ry-iy, wz = rz-iz;
    return (1-wx)*(1-wy)*(1-wz)*(*this)(ix,iy,iz)
      + wx*(1-wy)*(1-wz)*(*this)(ixp1,iy,iz)
      + (1-wx)*wy*(1-wz)*(*this)(ix,iyp1,iz)
      + (1-wx)*(1-wy)*wz*(*this)(ix,iy,izp1)
      + wx*(1-wy)*wz*(*this)(ixp1,iy,izp1)
      + (1-wx)*wy*wz*(*this)(ix,iyp1,izp1)
      + wx*wy*(1-wz)*(*this)(ixp1,iyp1,iz)
      + wx*wy*wz*(*this)(ixp1,iyp1,izp1);
  }
private:
  Lattice Lat;
  int Nx, Ny, Nz, NyNz, NxNyNz;
};
