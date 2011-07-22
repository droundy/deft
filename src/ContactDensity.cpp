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

#include "ContactDensity.h"

static Functional VectorThirdTerm(double radius) {
  Functional R(radius, "R");
  Functional n2 = ShellConvolve(radius);
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  return n2*(sqr(n2) - 3*(sqr(n2x) + sqr(n2y) + sqr(n2z)));
}

static Functional TensorThirdTerm(double radius) {
  Functional R(radius, "R");
  Functional n2 = ShellConvolve(radius);
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  Functional nTxx = xxShellConvolve(radius);
  Functional nTyy = yyShellConvolve(radius);
  Functional nTzz = zzShellConvolve(radius);
  Functional nTxy = xyShellConvolve(radius);
  Functional nTyz = yzShellConvolve(radius);
  Functional nTzx = zxShellConvolve(radius);
  Functional nTxz = nTzx, nTyx = nTxy, nTzy = nTyz;
  /*
  Functional trace_nT3 =
    nTxx*nTxx*nTxx + nTxx*nTxy*nTyx + nTxx*nTxz*nTzx + // starting with nTxx
    nTxy*nTyx*nTxx + nTxy*nTyy*nTyx + nTxy*nTyz*nTzx + // starting with nTxy
    nTxz*nTzx*nTxx + nTxz*nTzy*nTyx + nTxz*nTzz*nTzx + // starting with nTxz

    nTyx*nTxx*nTxy + nTyx*nTxy*nTyy + nTyx*nTxz*nTzy + // starting with nTyx
    nTyy*nTyx*nTxy + nTyy*nTyy*nTyy + nTyy*nTyz*nTzy + // starting with nTyy
    nTyz*nTzx*nTxy + nTyz*nTzy*nTyy + nTyz*nTzz*nTzy + // starting with nTyz

    nTzx*nTxx*nTxz + nTzx*nTxy*nTyz + nTzx*nTxz*nTzz + // starting with nTzx
    nTzy*nTyx*nTxz + nTzy*nTyy*nTyz + nTzy*nTyz*nTzz + // starting with nTzy
    nTzz*nTzx*nTxz + nTzz*nTzy*nTyz + nTzz*nTzz*nTzz; // starting with nTzz
  */
  Functional trace_nT3 =
    6*nTxy*nTyz*nTzx +
    nTxx*(  sqr(nTxx) + 3*sqr(nTxy) + 3*sqr(nTzx)) +
    nTyy*(3*sqr(nTxy) +   sqr(nTyy) + 3*sqr(nTyz)) +
    nTzz*(3*sqr(nTzx) + 3*sqr(nTzy) +   sqr(nTzz));
  return (n2*sqr(n2) - 3*n2*(sqr(n2x) + sqr(n2y) + sqr(n2z))
          +
          9*(sqr(n2x)*nTxx + sqr(n2y)*nTyy + sqr(n2z)*nTzz
             + 2*(n2x*n2y*nTxy + n2y*n2z*nTyz + n2z*n2x*nTzx)
             - 0.5*trace_nT3));
}

static Functional n3(double radius) {
  return StepConvolve(radius);
}

static Functional n2(double radius) {
  return ShellConvolve(radius);
}

static Functional n1(double radius) {
  Functional R(radius, "R");
  return ShellConvolve(radius)/(4*M_PI*R);
}

static Functional n0(double radius) {
  Functional R(radius, "R");
  return ShellConvolve(radius)/(4*M_PI*sqr(R));
}

Functional dWBNT_dn0(double radius) {
  return -log(1-n3(radius));
}

Functional dWBNT_dn1(double radius) {
  return n2(radius)/(1-n3(radius));
}

Functional dWBNT_dn2(double radius) {
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  return n1(radius)/(1-n3(radius)) +
    (sqr(n2(radius)) - (sqr(n2x) + sqr(n2y) + sqr(n2z)))*
    (n3(radius) + sqr(1-n3(radius))*log(1-n3(radius)))/
    (12*M_PI*sqr(n3(radius))*sqr(1-n3(radius)));
}

Functional dWBNT_dn3(double radius) {
  Functional R(radius, "R");
  Functional vtt = VectorThirdTerm(radius);
  Functional n3_ = n3(radius);
  Functional omn3 = 1 - n3_;
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  Functional nV22 = sqr(n2x) + sqr(n2y) + sqr(n2z);
  Functional n22mnV22 = sqr(n2(radius)) - nV22;
  return n0(radius)/omn3 +
    n22mnV22/sqr(omn3)/(4*M_PI*R) +
    (1/(36*M_PI))*vtt*(Functional(2)/(n3_*Pow(3)(omn3))
         - Functional(1)/(sqr(n3_)*sqr(omn3))
         - Functional(1)/(sqr(n3_)*omn3)
         - 2*log(omn3)/Pow(3)(n3_));
  return n0(radius)/omn3 +
    R*n22mnV22/sqr(omn3) +
    // vtt*(Functional(1)/(n3_*Pow(3)(omn3))
    //      - Functional(1)/(sqr(n3_)*sqr(omn3))
    //      - Functional(1)/(sqr(n3_)*omn3)
    //      - 2*log(omn3)/Pow(3)(n3_));
    vtt*(Functional(1)/sqr(omn3) + 2*log(omn3)/(35*M_PI*sqr(n3_)*sqr(omn3)));
  Functional ttt = TensorThirdTerm(radius);
}

Functional dAdR_simplest(double radius) {
  Functional R(radius, "R");
  return dWBNT_dn3(radius)*n2(radius) +
    dWBNT_dn2(radius)*(Functional(2)/R)*n2(radius) +
    dWBNT_dn1(radius)*n0(radius);
}

Functional dAdR_sphere_over_n(double radius) {
  Functional R(radius, "R");
  Functional two_over_R = Functional(2)/R;
  Functional one_over_4piRsqr = Functional(1)/(4*M_PI*sqr(R));
  return ShellConvolve(radius)(dWBNT_dn3(radius)
                               + two_over_R*dWBNT_dn2(radius)
                               + one_over_4piRsqr*dWBNT_dn1(radius));
}

Functional ContactDensitySimplest(double radius) {
  return Functional(0.25)*dAdR_simplest(radius)/n2(radius);
}

Functional ContactDensitySphere(double radius) {
  Functional R(radius, "R");
  return (Functional(1)/(4*(4*M_PI*sqr(R))))*dAdR_sphere_over_n(radius);
}

Functional FuWuContactDensity(double radius) {
  Functional R(radius, "R");
  Functional n2 = ShellConvolve(radius);
  Functional n0 = n2/(4*M_PI*sqr(R));
  Functional n3 = StepConvolve(radius);
  // The following is from equation 13 of Fu and Wu 2005, which I have
  // translated in terms of n2.  zeta3 is a version of the packing
  // fraction (usually called eta in our code) that is computed using
  // the shell convolution, so it is using the weighted density that
  // is more direcly relevant to the association free energy.
  Functional zeta3 = (R/Functional(3))*n2;

  // This gHS (called simply gHS) is the gHS that is used in Fu and
  // Wu's 2005 paper, in equation 13.  It seems ideal, since it
  // includes spatial dependence and is published and tested in
  // various ways.

  Functional ghs = gHS(zeta3, radius);
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  Functional zeta = 1 - (sqr(n2x) + sqr(n2y) + sqr(n2z))/sqr(n2);
  return n0*zeta*ghs;
}
