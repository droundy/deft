// -*- mode: C++; -*-

#pragma once

#include "lattice.h"
#include <stdio.h>
#include <eigen2/Eigen/Geometry>

static const double default_eps_size = 1000.0;

template<typename Scalar>
struct any_op EIGEN_EMPTY_STRUCT {
  EIGEN_STRONG_INLINE any_op(Lattice l, int nx, int ny, int nz,
                             Scalar (*func)(Cartesian)) : lat(l), fun(func) {
    Nx = nx; Ny = ny; Nz = nz;
    NyNz = Ny*Nz; NxNyNz = Nx*NyNz;
    dx = 1.0/Nx; dy = 1.0/Ny; dz = 1.0/Nz;
  }
  EIGEN_STRONG_INLINE const Scalar operator() (int row, int col) const {
    int n = row + col;
    const int z = n % Nz;
    n = (n-z)/Nz;
    const int y = n % Ny;
    const int x = (n-y)/Ny;
    const Relative rvec(x*dx,y*dy,z*dz);
    return fun(lat.wignerSeitz(lat.toCartesian(rvec)));
  }
  Lattice lat;
  int Nx, Ny, Nz, NyNz, NxNyNz;
  double dx, dy, dz;
  Scalar (*fun)(Cartesian);
};
namespace Eigen {
  template<typename Scalar>
  struct ei_functor_traits<any_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
} // namespace Eigen

class GridDescription {
public:
  explicit GridDescription(Lattice lat, int nx, int ny, int nz);
  GridDescription(const GridDescription &x);

  Lattice Lat, fineLat;
  int Nx, Ny, Nz, NyNz, NxNyNz;
  double dx, dy, dz;
};

class Grid : public VectorXd {
public:
  explicit Grid(Lattice lat, int nx, int ny, int nz);
  Grid(const Grid &x);

  // We need to define this for our object to work
  typedef Eigen::VectorXd Base;
  template<typename OtherDerived>
  Grid &operator=(const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  double operator()(int x, int y, int z) const {
    return (*this)[x*gd.NyNz + y*gd.Nz + z];
  }
  double &operator()(int x, int y, int z) {
    return (*this)[x*gd.NyNz + y*gd.Nz + z];
  }
  double operator()(const Cartesian &r) const {
    return (*this)(gd.Lat.toRelative(r));
  }
  double operator()(const Relative &r) const;
  void Set(double f(Cartesian));
  void epsSlice(const char *fname, Cartesian xmax, Cartesian ymax,
                Cartesian corner, int resolution) const;
  void epsNativeSlice(const char *fname,
                      Cartesian xmax, Cartesian ymax, Cartesian corner) const;
  Eigen::CwiseNullaryOp<any_op<double>, VectorXd> r2() const {
    return NullaryExpr(gd.NxNyNz, 1, r2_op);
  }
  Eigen::CwiseNullaryOp<any_op<double>, VectorXd> x() const {
    return NullaryExpr(gd.NxNyNz, 1, x_op);
  }
  Eigen::CwiseNullaryOp<any_op<double>, VectorXd> y() const {
    return NullaryExpr(gd.NxNyNz, 1, y_op);
  }
  Eigen::CwiseNullaryOp<any_op<double>, VectorXd> z() const {
    return NullaryExpr(gd.NxNyNz, 1, z_op);
  }
private:
  GridDescription gd;
  any_op<double> r2_op, x_op, y_op, z_op;
};

class ReciprocalGrid : public VectorXd {
public:
  explicit ReciprocalGrid(Lattice lat, int nx, int ny, int nz);
  ReciprocalGrid(const ReciprocalGrid &x);

  // We need to define this for our object to work
  typedef Eigen::VectorXd Base;
  template<typename OtherDerived>
  ReciprocalGrid &operator=(const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  double operator()(int x, int y, int z) const {
    return (*this)[x*gd.NyNz + y*gd.Nz + z];
  }
  double &operator()(int x, int y, int z) {
    return (*this)[x*gd.NyNz + y*gd.Nz + z];
  }
  double operator()(const Reciprocal &r) const {
    return (*this)(gd.Lat.toRelativeReciprocal(r));
  }
  double operator()(const RelativeReciprocal &) const {
    return 0;
    //return interpolate(r(0), r(1), r(2));
  }
  //void Set(double f(Cartesian));
  //void epsSlice(const char *fname, Cartesian xmax, Cartesian ymax,
  //              Cartesian corner, int resolution) const;
  //void epsNativeSlice(const char *fname,
  //                    Cartesian xmax, Cartesian ymax, Cartesian corner) const;
  Eigen::CwiseNullaryOp<any_op<double>, VectorXd> k2() const {
    return NullaryExpr(gd.NxNyNz, 1, r2_op);
  }
  Eigen::CwiseNullaryOp<any_op<double>, VectorXd> kx() const {
    return NullaryExpr(gd.NxNyNz, 1, x_op);
  }
  Eigen::CwiseNullaryOp<any_op<double>, VectorXd> ky() const {
    return NullaryExpr(gd.NxNyNz, 1, y_op);
  }
  Eigen::CwiseNullaryOp<any_op<double>, VectorXd> kz() const {
    return NullaryExpr(gd.NxNyNz, 1, z_op);
  }
private:
  GridDescription gd;
  any_op<double> r2_op, x_op, y_op, z_op;
};

/*
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
*/
