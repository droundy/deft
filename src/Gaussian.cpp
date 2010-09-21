// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Deft Authors
//
// Deft is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with deft; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include "Functionals.h"
#include "ReciprocalGrid.h"
#include <stdio.h>

class GaussianType : public FunctionalInterface {
public:
  GaussianType(double w) {
    width = w;
    kfac = -0.5*width*width; // FIXME: get width right in k space!
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    Grid out(gd, data);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= (kfac*g2(gd)).cwise().exp();
    return recip.ifft();
  }
  double transform(double n) const {
    return n;
  }
  double derive(double) const {
    return 1;
  }
  Functional grad(const Functional &ingrad, bool) const {
    return Gaussian(width)(ingrad);
  }
  void grad(const GridDescription &gd, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= (kfac*g2(gd)).cwise().exp();
    out = recip.ifft();
    *outgrad += out;

    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
  Expression printme(const Expression &x) const {
    return funexpr("Gaussian", Expression("width"))(Expression("gd"), x);
  }
private:
  double width, kfac;
};

Functional Gaussian(double width) {
  return Functional(new GaussianType(width));
}

static double myR, mydr;

Functional StepConvolve(double R) {
  return Functional(function_for_convolve<step_op<complex> >, R, true);
}

Functional ShellConvolve(double R) {
  return Functional(function_for_convolve<shell_op<complex> >, R, true);
}

static complex xdelta(Reciprocal kvec) {
  double k = kvec.norm();
  double kR = k*myR;
  if (kR > 1e-3) {
    return complex(0,exp(-spreading*k*k*mydr*mydr)*(4*M_PI)*kvec[0]/(k*k)*(myR*cos(kR) - sin(kR)/k));
  } else {
    const double kR2 = kR*kR;
    // The following is a simple power series expansion to the above
    // function, to handle the case as k approaches zero with greater
    // accuracy (and efficiency).  I evaluate the smaller elements
    // first in the hope of reducing roundoff error (but this is not
    // yet tested).
    return complex(0,(4*M_PI)*myR*myR*myR*kvec[0]*
                   (kR2*kR2*(1.0/5040-1.0/720)+kR2*(1.0/24-1.0/120)-1.0/3));
  }
}

static complex ydelta(Reciprocal kvec) {
  double k = kvec.norm();
  double kR = k*myR;
  if (kR > 1e-3) {
    return complex(0,exp(-spreading*k*k*mydr*mydr)*(4*M_PI)*kvec[1]/(k*k)*(myR*cos(kR) - sin(kR)/k));
  } else {
    const double kR2 = kR*kR;
    // The following is a simple power series expansion to the above
    // function, to handle the case as k approaches zero with greater
    // accuracy (and efficiency).  I evaluate the smaller elements
    // first in the hope of reducing roundoff error (but this is not
    // yet tested).
    return complex(0,(4*M_PI)*myR*myR*myR*kvec[1]*
                   (kR2*kR2*(1.0/5040-1.0/720)+kR2*(1.0/24-1.0/120)-1.0/3));
  }
}

static complex zdelta(Reciprocal kvec) {
  double k = kvec.norm();
  double kR = k*myR;
  if (kR > 1e-3) {
    return complex(0,exp(-spreading*k*k*mydr*mydr)*(4*M_PI)*kvec[2]/(k*k)*(myR*cos(kR) - sin(kR)/k));
  } else {
    const double kR2 = kR*kR;
    // The following is a simple power series expansion to the above
    // function, to handle the case as k approaches zero with greater
    // accuracy (and efficiency).  I evaluate the smaller elements
    // first in the hope of reducing roundoff error (but this is not
    // yet tested).
    return complex(0,(4*M_PI)*myR*myR*myR*kvec[2]*
                   (kR2*kR2*(1.0/5040-1.0/720)+kR2*(1.0/24-1.0/120)-1.0/3));
  }
}

class VShellConvolveType : public FunctionalInterface {
public:
  VShellConvolveType(double radius, int dir) : R(radius), direction(dir) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    ReciprocalGrid recip(gd);
    {
      const Grid out(gd, data);
      recip = out.fft();
      // We want to free out immediately to save memory!
    }
    myR = R;
    mydr = pow(gd.fineLat.volume(), 1.0/3);
    switch (direction) {
    case 0: recip.MultiplyBy(xdelta);
      break;
    case 1: recip.MultiplyBy(ydelta);
      break;
    case 2: recip.MultiplyBy(zdelta);
      break;
    }
    return recip.ifft();
  }
  double transform(double) const {
    return 0;
  }
  double derive(double) const {
    return 0;
  }
  Functional grad(const Functional &ingrad, bool) const {
    return Functional(new VShellConvolveType(R, direction))((-1)*ingrad);
  }
  void grad(const GridDescription &gd, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    myR = R;
    mydr = pow(gd.fineLat.volume(), 1.0/3);
    switch (direction) {
    case 0: recip.MultiplyBy(xdelta);
      break;
    case 1: recip.MultiplyBy(ydelta);
      break;
    case 2: recip.MultiplyBy(zdelta);
      break;
    }
    out = recip.ifft();
    *outgrad -= out;

    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
  Expression printme(const Expression &x) const {
    switch (direction) {
    case 0: return funexpr("xShellConvolve", Expression("R"))(Expression("gd"), x);
    case 1: return funexpr("yShellConvolve", Expression("R"))(Expression("gd"), x);
    default:
      return funexpr("zShellConvolve", Expression("R"))(Expression("gd"), x);
    }
  }
private:
  double R;
  int direction;
};

Functional xShellConvolve(double R) {
  //return Functional(function_for_convolve<xshell_op<complex> >, R, false);
  return Functional(new VShellConvolveType(R, 0));
}
Functional yShellConvolve(double R) {
  //return Functional(function_for_convolve<yshell_op<complex> >, R, false);
  return Functional(new VShellConvolveType(R, 1));
}
Functional zShellConvolve(double R) {
  //return Functional(function_for_convolve<zshell_op<complex> >, R, false);
  return Functional(new VShellConvolveType(R, 2));
}

Functional xyShellConvolve(double R) {
  return Functional(function_for_convolve<xyshell_op<complex> >, R, true);
}
Functional yzShellConvolve(double R) {
  return Functional(function_for_convolve<yzshell_op<complex> >, R, true);
}
Functional zxShellConvolve(double R) {
  return Functional(function_for_convolve<zxshell_op<complex> >, R, true);
}

Functional xxShellConvolve(double R) {
  return Functional(function_for_convolve<xxshell_op<complex> >, R, true);
}
Functional yyShellConvolve(double R) {
  return Functional(function_for_convolve<yyshell_op<complex> >, R, true);
}
Functional zzShellConvolve(double R) {
  return Functional(function_for_convolve<zzshell_op<complex> >, R, true);
}
