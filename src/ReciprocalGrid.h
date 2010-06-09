// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include <stdio.h>
#include <eigen2/Eigen/Geometry>

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


class ReciprocalGrid : public VectorXd {
public:
  explicit ReciprocalGrid(const GridDescription &);
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
  double operator()(const RelativeReciprocal &) const;
  Eigen::CwiseNullaryOp<any_rop<double>, VectorXd> k2() const {
    return NullaryExpr(gd.NxNyNz, 1, r2_op);
  }
  Eigen::CwiseNullaryOp<any_rop<double>, VectorXd> kx() const {
    return NullaryExpr(gd.NxNyNz, 1, x_op);
  }
  Eigen::CwiseNullaryOp<any_rop<double>, VectorXd> ky() const {
    return NullaryExpr(gd.NxNyNz, 1, y_op);
  }
  Eigen::CwiseNullaryOp<any_rop<double>, VectorXd> kz() const {
    return NullaryExpr(gd.NxNyNz, 1, z_op);
  }
private:
  GridDescription gd;
  any_rop<double> r2_op, x_op, y_op, z_op;
};
