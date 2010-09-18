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
  template<typename Scalar>
  struct ei_functor_traits<step_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<shell_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<xshell_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<yshell_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<zshell_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<xyshell_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<yzshell_op<Scalar> > {
    enum {
      Cost = NumTraits<Scalar>::AddCost,
      PacketAccess = false,
      IsRepeatable = true
    };
  };
  template<typename Scalar>
  struct ei_functor_traits<zxshell_op<Scalar> > {
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

// step
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<step_op<complex>, VectorXcd> step(const GridDescription &gd, double R) {
  return Eigen::CwiseNullaryOp<step_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, step_op<complex>(gd,R));
}
EIGEN_STRONG_INLINE step_op<complex> stepper(const GridDescription &gd, double R) {
  return step_op<complex>(gd,R);
}

// shell
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<shell_op<complex>, VectorXcd> shell(const GridDescription &gd, double R) {
  return Eigen::CwiseNullaryOp<shell_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, shell_op<complex>(gd,R));
}
EIGEN_STRONG_INLINE shell_op<complex> sheller(const GridDescription &gd, double R) {
  return shell_op<complex>(gd,R);
}

// x shell
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<xshell_op<complex>, VectorXcd> xshell(const GridDescription &gd, double R) {
  return Eigen::CwiseNullaryOp<xshell_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, xshell_op<complex>(gd,R));
}
EIGEN_STRONG_INLINE xshell_op<complex> xsheller(const GridDescription &gd, double R) {
  return xshell_op<complex>(gd,R);
}

// y shell
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<yshell_op<complex>, VectorXcd> yshell(const GridDescription &gd, double R) {
  return Eigen::CwiseNullaryOp<yshell_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, yshell_op<complex>(gd,R));
}
EIGEN_STRONG_INLINE yshell_op<complex> ysheller(const GridDescription &gd, double R) {
  return yshell_op<complex>(gd,R);
}

// z shell
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<zshell_op<complex>, VectorXcd> zshell(const GridDescription &gd, double R) {
  return Eigen::CwiseNullaryOp<zshell_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, zshell_op<complex>(gd,R));
}
EIGEN_STRONG_INLINE zshell_op<complex> zsheller(const GridDescription &gd, double R) {
  return zshell_op<complex>(gd,R);
}

// xy shell
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<xyshell_op<complex>, VectorXcd> xyshell(const GridDescription &gd, double R) {
  return Eigen::CwiseNullaryOp<xyshell_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, xyshell_op<complex>(gd,R));
}
EIGEN_STRONG_INLINE xyshell_op<complex> xysheller(const GridDescription &gd, double R) {
  return xyshell_op<complex>(gd,R);
}

// yz shell
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<yzshell_op<complex>, VectorXcd> yzshell(const GridDescription &gd, double R) {
  return Eigen::CwiseNullaryOp<yzshell_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, yzshell_op<complex>(gd,R));
}
EIGEN_STRONG_INLINE yzshell_op<complex> yzsheller(const GridDescription &gd, double R) {
  return yzshell_op<complex>(gd,R);
}

// zx shell
EIGEN_STRONG_INLINE Eigen::CwiseNullaryOp<zxshell_op<complex>, VectorXcd> zxshell(const GridDescription &gd, double R) {
  return Eigen::CwiseNullaryOp<zxshell_op<complex>, VectorXcd>(gd.NxNyNzOver2, 1, zxshell_op<complex>(gd,R));
}
EIGEN_STRONG_INLINE zxshell_op<complex> zxsheller(const GridDescription &gd, double R) {
  return zxshell_op<complex>(gd,R);
}
