// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include <stdio.h>
#include <eigen2/Eigen/Geometry>

typedef std::complex<double> complex;

template<typename Scalar>
struct any_rop EIGEN_EMPTY_STRUCT {
  EIGEN_STRONG_INLINE any_rop(const GridDescription &gdin,
                              Scalar (*func)(Reciprocal))
    : gd(gdin), fun(func) {}
  EIGEN_STRONG_INLINE const Scalar operator() (int row, int col) const {
    int n = row + col;
    const int z = n % gd.Nz;
    n = (n-z)/gd.Nz;
    const int y = n % gd.Ny;
    const int x = (n-y)/gd.Ny;
    const RelativeReciprocal rvec(x*gd.dx,y*gd.dy,z*gd.dz);
    return fun(gd.Lat.brillouinZone(gd.Lat.toReciprocal(rvec)));
  }
  GridDescription gd;
  Scalar (*fun)(Reciprocal);
};
namespace Eigen {
  template<typename Scalar>
  struct ei_functor_traits<any_rop<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
} // namespace Eigen

class Grid;

class ReciprocalGrid : public VectorXcd {
public:
  explicit ReciprocalGrid(const GridDescription &);
  ReciprocalGrid(const ReciprocalGrid &x);

  // We need to define this for our object to work
  typedef Eigen::VectorXcd Base;
  template<typename OtherDerived>
  ReciprocalGrid &operator=(const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  Grid ifft() const;
  complex operator()(int x, int y, int z) const {
    if (z > gd.NzOver2) z = gd.Nz - z;
    return (*this)[x*gd.NyNz + y*gd.Nz + z];
  }
  complex &operator()(int x, int y, int z) {
    if (z > gd.NzOver2) z = gd.Nz - z;
    return (*this)[x*gd.NyNz + y*gd.Nz + z];
  }
  complex operator()(const Reciprocal &r) const {
    return (*this)(gd.Lat.toRelativeReciprocal(r));
  }
  complex operator()(const RelativeReciprocal &) const;
  Eigen::CwiseNullaryOp<any_rop<double>, VectorXcd> g2() const {
    return NullaryExpr(gd.NxNyNz, 1, g2_op);
  }
  Eigen::CwiseNullaryOp<any_rop<double>, VectorXcd> gx() const {
    return NullaryExpr(gd.NxNyNz, 1, gx_op);
  }
  Eigen::CwiseNullaryOp<any_rop<double>, VectorXcd> gy() const {
    return NullaryExpr(gd.NxNyNz, 1, gy_op);
  }
  Eigen::CwiseNullaryOp<any_rop<double>, VectorXcd> gz() const {
    return NullaryExpr(gd.NxNyNz, 1, gz_op);
  }
private:
  GridDescription gd;
  any_rop<double> g2_op, gx_op, gy_op, gz_op;
};
