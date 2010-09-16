// -*- mode: C++; -*-

// This file contains "nullary" operators in reciprocal space.  These
// take very little memory (they just involve an additional copy of
// the GridDescription, since I'm lazy, but that's small), and allow
// us to treat things like "g^2" and "g_x" (almost) as VectorXcd's
// that take up no space, although they do take some time to work out
// what the position is for each grid point.  They are complex (even
// though the vectors themselves are real), because fields in fourier
// space are generally complex, so it makes it easier to mix and
// match.  If it turns out that eigen is smart about double-complex
// operations (which can be done more efficiently than converting the
// double to a complex), then this should be rewritten for greater
// speed.

#pragma once

#include "GridDescription.h"

template<typename Scalar>
struct base_rop : private GridDescription {
  EIGEN_STRONG_INLINE base_rop(const GridDescription &gdin) : GridDescription(gdin) {}
  EIGEN_STRONG_INLINE const Scalar operator() (int row, int col) const {
    int n = row + col;
    const int z = n % NzOver2;
    n = (n-z)/NzOver2;
    const int y = n % Ny;
    const int x = (n-y)/Ny;
    const RelativeReciprocal rvec((x>Nx/2) ? x - Nx : x,
                                  (y>Ny/2) ? y - Ny : y,
                                  z);
    // FIXME: it seems that brillouinZone is broken...  :(
    //return func(fineLat.brillouinZone(Lat.toReciprocal(rvec)));
    return func(Lat.toReciprocal(rvec));
  }
  virtual Scalar func(Reciprocal) const = 0;
};

template<typename Scalar>
struct g2_op : public base_rop<Scalar> {
  g2_op(const GridDescription &gd) : base_rop<Scalar>(gd) {}
  Scalar func(Reciprocal g) const {
    return g.squaredNorm();
  }
};

template<typename Scalar>
struct gx_op : public base_rop<Scalar> {
  gx_op(const GridDescription &gd) : base_rop<Scalar>(gd) {}
  Scalar func(Reciprocal g) const {
    return g(0);
  }
};

template<typename Scalar>
struct gy_op : public base_rop<Scalar> {
  gy_op(const GridDescription &gd) : base_rop<Scalar>(gd) {}
  Scalar func(Reciprocal g) const {
    return g(1);
  }
};

template<typename Scalar>
struct gz_op : public base_rop<Scalar> {
  gz_op(const GridDescription &gd) : base_rop<Scalar>(gd) {}
  Scalar func(Reciprocal g) const {
    return g(2);
  }
};

namespace Eigen {
  template<typename Scalar>
  struct ei_functor_traits<g2_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<gx_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<gy_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<gz_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true 
    };
  };
} // namespace Eigen

EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<g2_op<complex>, VectorXcd> g2(const GridDescription &gd) {
  return Eigen::CwiseNullaryOp<g2_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, g2_op<complex>(gd));
}

EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<gx_op<complex>, VectorXcd> gx(const GridDescription &gd) {
  return Eigen::CwiseNullaryOp<gx_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, gx_op<complex>(gd));
}

EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<gy_op<complex>, VectorXcd> gy(const GridDescription &gd) {
  return Eigen::CwiseNullaryOp<gy_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, gy_op<complex>(gd));
}

EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<gz_op<complex>, VectorXcd> gz(const GridDescription &gd) {
  return Eigen::CwiseNullaryOp<gz_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, gz_op<complex>(gd));
}
