// -*- mode: C++; -*-

// This file contains "nullary" operators in real space.  These take
// very little memory (they just involve an additional copy of the
// GridDescription, since I'm lazy, but that's small), and allow us to
// treat things like "r^2" and "r_x" (almost) as VectorXd's that take
// up no space, although they do take some time to work out what the
// position is for each grid point.

#pragma once

#include "GridDescription.h"

template<typename Scalar>
struct base_op : private GridDescription {
  EIGEN_STRONG_INLINE base_op(const GridDescription &gdin) : GridDescription(gdin) {}
  EIGEN_STRONG_INLINE const Scalar operator() (int row, int col) const {
    int n = row + col;
    const int z = n % Nz;
    n = (n-z)/Nz;
    const int y = n % Ny;
    const int x = (n-y)/Ny;
    const Relative rvec(x*dx,y*dy,z*dz);
    return func(Lat.wignerSeitz(Lat.toCartesian(rvec)));
  }
  virtual Scalar func(Cartesian) const = 0;
};

template<typename Scalar>
struct r2_op : public base_op<Scalar> {
  r2_op(const GridDescription &gd) : base_op<Scalar>(gd) {}
  Scalar func(Cartesian g) const {
    return g.squaredNorm();
  }
};

template<typename Scalar>
struct rx_op : public base_op<Scalar> {
  rx_op(const GridDescription &gd) : base_op<Scalar>(gd) {}
  Scalar func(Cartesian g) const {
    return g(0);
  }
};

template<typename Scalar>
struct ry_op : public base_op<Scalar> {
  ry_op(const GridDescription &gd) : base_op<Scalar>(gd) {}
  Scalar func(Cartesian g) const {
    return g(1);
  }
};

template<typename Scalar>
struct rz_op : public base_op<Scalar> {
  rz_op(const GridDescription &gd) : base_op<Scalar>(gd) {}
  Scalar func(Cartesian g) const {
    return g(2);
  }
};

namespace Eigen {
  template<typename Scalar>
  struct ei_functor_traits<r2_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<rx_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<ry_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<rz_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
} // namespace Eigen

EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<r2_op<double>, VectorXd> r2(const GridDescription &gd) {
  return Eigen::CwiseNullaryOp<r2_op<double>, VectorXd>(gd.NxNyNz, 1, r2_op<double>(gd));
}

EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<rx_op<double>, VectorXd> rx(const GridDescription &gd) {
  return Eigen::CwiseNullaryOp<rx_op<double>, VectorXd>(gd.NxNyNz, 1, rx_op<double>(gd));
}

EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<ry_op<double>, VectorXd> ry(const GridDescription &gd) {
  return Eigen::CwiseNullaryOp<ry_op<double>, VectorXd>(gd.NxNyNz, 1, ry_op<double>(gd));
}

EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<rz_op<double>, VectorXd> rz(const GridDescription &gd) {
  return Eigen::CwiseNullaryOp<rz_op<double>, VectorXd>(gd.NxNyNz, 1, rz_op<double>(gd));
}
