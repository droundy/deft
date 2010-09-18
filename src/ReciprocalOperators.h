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
    //const RelativeReciprocal rvec((x+Nx/2) % Nx - Nx/2,
    //                              (y+Ny/2) % Ny - Ny/2,
    //                              z);
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

static const double spreading = 3.0;

template<typename Scalar>
struct step_op : public base_rop<Scalar> {
  step_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    R3 = R*R*R;
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double kR = k*R;
    double kdr = k*dr;
    if (kR > 1e-3) {
      return exp(-spreading*kdr*kdr)*(4*M_PI)*(sin(kR) - kR*cos(kR))/(k*k*k);
    } else {
      const double kR2 = kR*kR;
      // The following is a simple power series expansion to the above
      // function, to handle the case as k approaches zero with greater
      // accuracy (and efficiency).  I evaluate the smaller elements
      // first in the hope of reducing roundoff error (but this is not
      // yet tested).
      return (4*M_PI/3)*R3*(kR2*kR2*kR2*(-1.0/15120) + kR2*kR2*(1.0/280) + kR2*(-1.0/10) + 1 );
    }
  }
  double R, R3, dr;
};

template<typename Scalar>
struct shell_op : public base_rop<Scalar> {
  shell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double kR = k*R;
    if (kR > 1e-3) {
      return exp(-spreading*k*k*dr*dr)*(4*M_PI)*R*sin(kR)/k;
    } else {
      const double kR2 = kR*kR;
      // The following is a simple power series expansion to the above
      // function, to handle the case as k approaches zero with greater
      // accuracy (and efficiency).  I evaluate the smaller elements
      // first in the hope of reducing roundoff error (but this is not
      // yet tested).
      return (4*M_PI)*(R*R)*(- kR2*kR2*kR2*(1.0/120/6/7) + kR2*kR2*(1.0/120) - kR2*(1.0/6) + 1);
  }
  }
  double R, dr;
};

template<typename Scalar>
struct xshell_op : public base_rop<Scalar> {
  xshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double kR = k*R;
    if (kR > 1e-3) {
      return complex(0,exp(-spreading*k*k*dr*dr)*(4*M_PI)*kvec[0]/(k*k)*(R*cos(kR) - sin(kR)/k));
    } else {
      const double kR2 = kR*kR;
      // The following is a simple power series expansion to the above
      // function, to handle the case as k approaches zero with greater
      // accuracy (and efficiency).  I evaluate the smaller elements
      // first in the hope of reducing roundoff error (but this is not
      // yet tested).
      return complex(0,(4*M_PI)*R*R*R*kvec[0]*
                     (kR2*kR2*(1.0/5040-1.0/720)+kR2*(1.0/24-1.0/120)-1.0/3));
    }
  }
  double R, dr;
};

template<typename Scalar>
struct yshell_op : public base_rop<Scalar> {
  yshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double kR = k*R;
    if (kR > 1e-3) {
      return complex(0,exp(-spreading*k*k*dr*dr)*(4*M_PI)*kvec[1]/(k*k)*(R*cos(kR) - sin(kR)/k));
    } else {
      const double kR2 = kR*kR;
      // The following is a simple power series expansion to the above
      // function, to handle the case as k approaches zero with greater
      // accuracy (and efficiency).  I evaluate the smaller elements
      // first in the hope of reducing roundoff error (but this is not
      // yet tested).
      return complex(0,(4*M_PI)*R*R*R*kvec[1]*
                     (kR2*kR2*(1.0/5040-1.0/720)+kR2*(1.0/24-1.0/120)-1.0/3));
    }
  }
  double R, dr;
};

template<typename Scalar>
struct zshell_op : public base_rop<Scalar> {
  zshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double kR = k*R;
    if (kR > 1e-3) {
      return complex(0,exp(-spreading*k*k*dr*dr)*(4*M_PI)*kvec[2]/(k*k)*(R*cos(kR) - sin(kR)/k));
    } else {
      const double kR2 = kR*kR;
      // The following is a simple power series expansion to the above
      // function, to handle the case as k approaches zero with greater
      // accuracy (and efficiency).  I evaluate the smaller elements
      // first in the hope of reducing roundoff error (but this is not
      // yet tested).
      return complex(0,(4*M_PI)*R*R*R*kvec[2]*
                     (kR2*kR2*(1.0/5040-1.0/720)+kR2*(1.0/24-1.0/120)-1.0/3));
    }
  }
  double R, dr;
};

template<typename Scalar>
struct xyshell_op : public base_rop<Scalar> {
  xyshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double kR = k*R;
    if (kR > 0)
      return -exp(-spreading*k*k*dr*dr)*(4*M_PI)*kvec[0]*kvec[1]/(k*k*k*k)*((3-kR*kR)*sin(kR)/kR - 3*cos(kR));
    else
      return 0;
  }
  double R, dr;
};

template<typename Scalar>
struct yzshell_op : public base_rop<Scalar> {
  yzshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double kR = k*R;
    if (kR > 0)
      return -exp(-spreading*k*k*dr*dr)*(4*M_PI)*kvec[1]*kvec[2]/(k*k*k*k)*((3-kR*kR)*sin(kR)/kR - 3*cos(kR));
    else
      return 0;
  }
  double R, dr;
};

template<typename Scalar>
struct zxshell_op : public base_rop<Scalar> {
  zxshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double kR = k*R;
    if (kR > 0)
      return -exp(-spreading*k*k*dr*dr)*(4*M_PI)*kvec[2]*kvec[0]/(k*k*k*k)*((3-kR*kR)*sin(kR)/kR - 3*cos(kR));
    else
      return 0;
  }
  double R, dr;
};

template<typename Scalar>
struct xxshell_op : public base_rop<Scalar> {
  xxshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double k2 = k*k;
    double kR = k*R;
    if (kR > 0) {
      double cos2 = kvec[0]*kvec[0]/k2;
      return exp(-spreading*k*k*dr*dr)*
        (4*M_PI*(R*R))*(sin(kR)/kR*(-1./3 + cos2 - 3*cos2/(kR*kR) + 1/(kR*kR)) +
                        cos(kR)/(kR*kR)*(-1 + 3*cos2));
    } else
      return 0;
  }
  double R, dr;
};

template<typename Scalar>
struct yyshell_op : public base_rop<Scalar> {
  yyshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double k2 = k*k;
    double kR = k*R;
    if (kR > 0) {
      double cos2 = kvec[1]*kvec[1]/k2;
      return exp(-spreading*k*k*dr*dr)*
        (4*M_PI*(R*R))*(sin(kR)/kR*(-1./3 + cos2 - 3*cos2/(kR*kR) + 1/(kR*kR)) +
                        cos(kR)/(kR*kR)*(-1 + 3*cos2));
    } else
      return 0;
  }
  double R, dr;
};

template<typename Scalar>
struct zzshell_op : public base_rop<Scalar> {
  zzshell_op(const GridDescription &gd, double r) : base_rop<Scalar>(gd), R(r) {
    dr = pow(gd.fineLat.volume(), 1.0/3);
  }
  Scalar func(Reciprocal kvec) const {
    double k = kvec.norm();
    double k2 = k*k;
    double kR = k*R;
    if (kR > 0) {
      double cos2 = kvec[2]*kvec[2]/k2;
      return exp(-spreading*k*k*dr*dr)*
        (4*M_PI*(R*R))*(sin(kR)/kR*(-1./3 + cos2 - 3*cos2/(kR*kR) + 1/(kR*kR)) +
                        cos(kR)/(kR*kR)*(-1 + 3*cos2));
    } else
      return 0;
  }
  double R, dr;
};

namespace Eigen {
#define ADD_ONE_ROP(t) template<typename Scalar> \
                       struct ei_functor_traits<t<Scalar> > { \
                         enum { \
                           Cost = NumTraits<Scalar>::AddCost, \
                           PacketAccess = false,              \
                           IsRepeatable = true                \
                         }; \
                       }
  ADD_ONE_ROP(g2_op);
  ADD_ONE_ROP(gx_op);
  ADD_ONE_ROP(gy_op);
  ADD_ONE_ROP(gz_op);
  ADD_ONE_ROP(step_op);
  ADD_ONE_ROP(shell_op);
  ADD_ONE_ROP(xshell_op);
  ADD_ONE_ROP(yshell_op);
  ADD_ONE_ROP(zshell_op);

  ADD_ONE_ROP(xyshell_op);
  ADD_ONE_ROP(yzshell_op);
  ADD_ONE_ROP(zxshell_op);

  ADD_ONE_ROP(xxshell_op);
  ADD_ONE_ROP(yyshell_op);
  ADD_ONE_ROP(zzshell_op);
} // namespace Eigen

#define MAKE_FUNCTION_ROP(t,f) \
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<t<complex>, VectorXcd> f(const GridDescription &gd) { \
  return Eigen::CwiseNullaryOp<t<complex>, VectorXcd>(gd.NxNyNzOver2, 1, t<complex>(gd)); \
}

MAKE_FUNCTION_ROP(g2_op,g2)

MAKE_FUNCTION_ROP(gx_op,gx)
MAKE_FUNCTION_ROP(gy_op,gy)
MAKE_FUNCTION_ROP(gz_op,gz)

#define MAKE_FUNCTION_WITH_DOUBLE_ROP(t,f) \
  EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<t<complex>, VectorXcd> f(const GridDescription &gd, double R) { \
    return Eigen::CwiseNullaryOp<t<complex>, VectorXcd>(gd.NxNyNzOver2, 1, t<complex>(gd, R)); \
  }

MAKE_FUNCTION_WITH_DOUBLE_ROP(step_op,step)
MAKE_FUNCTION_WITH_DOUBLE_ROP(shell_op,shell)

MAKE_FUNCTION_WITH_DOUBLE_ROP(xshell_op,xshell)
MAKE_FUNCTION_WITH_DOUBLE_ROP(yshell_op,yshell)
MAKE_FUNCTION_WITH_DOUBLE_ROP(zshell_op,zshell)

MAKE_FUNCTION_WITH_DOUBLE_ROP(xyshell_op,xyshell)
MAKE_FUNCTION_WITH_DOUBLE_ROP(yzshell_op,yzshell)
MAKE_FUNCTION_WITH_DOUBLE_ROP(zxshell_op,zxshell)

MAKE_FUNCTION_WITH_DOUBLE_ROP(xxshell_op,xxshell)
MAKE_FUNCTION_WITH_DOUBLE_ROP(yyshell_op,yyshell)
MAKE_FUNCTION_WITH_DOUBLE_ROP(zzshell_op,zzshell)

template<typename T>
EIGEN_STRONG_INLINE T function_for_convolve(const GridDescription &gd, double R) {
  return T(gd,R);
}
