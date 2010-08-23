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

class GaussianType : public FieldFunctionalInterface {
public:
  GaussianType(double width) {
    kfac = -0.5*width*width; // FIXME: get width right in k space!
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    Grid out(gd, data);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= (kfac*recip.g2()).cwise().exp();
    return recip.ifft();
  }
  double transform(double n) const {
    return n;
  }

  void grad(const GridDescription &gd, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= (kfac*recip.g2()).cwise().exp();
    out = recip.ifft();
    *outgrad += out;

    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
private:
  double kfac;
};

FieldFunctional Gaussian(double width) {
  return FieldFunctional(new GaussianType(width));
}

static double myR;
static double step(Reciprocal kvec) {
  double k = kvec.norm();
  double kR = k*myR;
  if (kR > 1e-3) {
    return (4*M_PI)*(sin(kR) - kR*cos(kR))/(k*k*k);
  } else {
    const double kR2 = kR*kR;
    // The following is a simple power series expansion to the above
    // function, to handle the case as k approaches zero with greater
    // accuracy (and efficiency).  I evaluate the smaller elements
    // first in the hope of reducing roundoff error (but this is not
    // yet tested).
    return (4*M_PI)*(myR*myR*myR)*(kR2*kR2*kR2*(-1.0/45360) + kR2*kR2*(1.0/840) + kR2*(-1.0/30) + (1.0/3) );
  }
}

class StepConvolveType : public FieldFunctionalInterface {
public:
  StepConvolveType(double radius) : R(radius) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    ReciprocalGrid recip(gd);
    {
      const Grid out(gd, data);
      recip = out.fft();
      // We want to free out immediately to save memory!
    }
    myR = R;
    recip.MultiplyBy(step);
    return recip.ifft();
  }
  double transform(double n) const {
    return n*(4*M_PI)*R*R*R;
  }

  void grad(const GridDescription &gd, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    myR = R;
    recip.MultiplyBy(step);
    out = recip.ifft();
    *outgrad += out;

    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
private:
  double R;
};

FieldFunctional StepConvolve(double R) {
  return FieldFunctional(new StepConvolveType(R));
}

static double delta(Reciprocal kvec) {
  double k = kvec.norm();
  double kR = k*myR;
  if (kR > 1e-3) {
    return (4*M_PI)*myR*sin(kR)/k;
  } else {
    const double kR2 = kR*kR;
    // The following is a simple power series expansion to the above
    // function, to handle the case as k approaches zero with greater
    // accuracy (and efficiency).  I evaluate the smaller elements
    // first in the hope of reducing roundoff error (but this is not
    // yet tested).
    return (4*M_PI)*(myR*myR)*(- kR2*kR2*kR2*(1.0/120/6/7) + kR2*kR2*(1.0/120) - kR2*(1.0/6) + 1);
  }
}

class ShellConvolveType : public FieldFunctionalInterface {
public:
  ShellConvolveType(double radius) : R(radius) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    ReciprocalGrid recip(gd);
    {
      const Grid out(gd, data);
      recip = out.fft();
      // We want to free out immediately to save memory!
    }
    myR = R;
    recip.MultiplyBy(delta);
    return recip.ifft();
  }
  double transform(double n) const {
    return n*(4*M_PI)*R*R*R;
  }

  void grad(const GridDescription &gd, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    myR = R;
    recip.MultiplyBy(delta);
    out = recip.ifft();
    *outgrad += out;

    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
private:
  double R;
};

FieldFunctional ShellConvolve(double R) {
  return FieldFunctional(new ShellConvolveType(R));
}

static complex xdelta(Reciprocal kvec) {
  double k = kvec.norm();
  double kR = k*myR;
  if (kR > 1e-3) {
    return complex(0,(4*M_PI)*myR*myR*myR*kvec[0]/(k*k)*(myR*cos(kR) - sin(kR)/k));
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
    return complex(0,(4*M_PI)*myR*myR*myR*kvec[1]/(k*k)*(myR*cos(kR) - sin(kR)/k));
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
    return complex(0,(4*M_PI)*myR*myR*myR*kvec[2]/(k*k)*(myR*cos(kR) - sin(kR)/k));
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

class VShellConvolveType : public FieldFunctionalInterface {
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

  void grad(const GridDescription &gd, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    myR = R;
    switch (direction) {
    case 0: recip.MultiplyBy(xdelta);
      break;
    case 1: recip.MultiplyBy(ydelta);
      break;
    case 2: recip.MultiplyBy(zdelta);
      break;
    }
    out = recip.ifft();
    *outgrad += out;

    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
private:
  double R;
  int direction;
};

FieldFunctional xShellConvolve(double R) {
  return FieldFunctional(new VShellConvolveType(R, 0));
}
FieldFunctional yShellConvolve(double R) {
  return FieldFunctional(new VShellConvolveType(R, 1));
}
FieldFunctional zShellConvolve(double R) {
  return FieldFunctional(new VShellConvolveType(R, 2));
}
