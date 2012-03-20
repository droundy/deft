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

Functional n0(double radius) {
  Functional R(radius, "R");
  return ShellConvolve(radius)/(4*M_PI*sqr(R));
}

static Functional n_A(double radius) {
  Functional R(2*radius, "R");
  return ShellConvolve(2*radius,Expression("2*R"))/(16*M_PI*sqr(R));
}

Functional dWB_dn0(double radius) {
  return -log(1-n3(radius));
}

Functional dWB_dn1(double radius) {
  return n2(radius)/(1-n3(radius));
}

Functional dWB_dn2(double radius) {
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  return n1(radius)/(1-n3(radius)) +
    (sqr(n2(radius)) - (sqr(n2x) + sqr(n2y) + sqr(n2z)))*
    (n3(radius) + sqr(1-n3(radius))*log(1-n3(radius)))/
    (12*M_PI*sqr(n3(radius))*sqr(1-n3(radius)));
}

Functional dWB_dn2v_over_n2v(double radius) {
  Functional R(radius, "R");
  return (Functional(1)/(4*M_PI*R))/(1-n3(radius)) -
    6*n2(radius)*(n3(radius) + sqr(1-n3(radius))*log(1-n3(radius))/(36*M_PI*sqr(n3(radius))*sqr(1-n3(radius))));
}

Functional dWB_dn1v_over_n2v(double radius) {
  return Functional(1)/(1-n3(radius));
}

Functional dWB_dn3(double radius) {
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



Functional dWBm2_dn0(double radius) {
  return -log(1-n3(radius));
}

Functional phi2(const Functional &n3) {
  return (6-3*n3+6*(1-n3)*log(1-n3)/n3)/sqr(n3);
}

Functional phi3(const Functional &n3) {
  return (6 - 9*n3 + 6*sqr(n3) + 6*sqr(1-n3)*log(1-n3)/n3)/(4*sqr(n3));
}

static Functional sqr_n3_times_phi2_dn3(const Functional &n3) {
  return ((Functional(-6)+4*n3)*log(1-n3)-6*n3+sqr(n3))*Functional(3)/sqr(n3);
}

static Functional n3_times_phi3_dn3(const Functional &n3) {
  return (Functional(-6) + 5*n3 + (-2*sqr(n3)+8*n3-Functional(6))*log(1-n3)/n3)*(Functional(3)/(4*sqr(n3)));
}

Functional dWBm2_dn1v_over_n2v(double radius) {
  return -(Functional(1) + (1.0/9.0)*sqr(n3(radius))*phi2(n3(radius))) / (1-n3(radius));
}

Functional dWBm2_dn1(double radius) {
  return -n2(radius)*dWBm2_dn1v_over_n2v(radius);
}

Functional dWBm2_dn2(double radius) {
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  return (Functional(1) + (1/9.0)*sqr(n3(radius))*phi2(n3(radius)))*n1(radius) / (1-n3(radius)) +
    (1 - (4.0/9.0)*n3(radius)*phi3(n3(radius)))*
    (sqr(n2(radius)) - (sqr(n2x) + sqr(n2y) + sqr(n2z))) / (8*M_PI*sqr(1-n3(radius)));
}

Functional dWBm2_dn2v_over_n2v(double radius) {
  Functional R(radius, "R");
  return -(Functional(1) + (1.0/9.0)*sqr(n3(radius))*phi2(n3(radius))) / (4*M_PI*R*(1-n3(radius))) 
    -(1 - (4.0/9.0)*n3(radius)*phi3(n3(radius)))*
    Functional(6)*n2(radius) / (24*M_PI*sqr(1-n3(radius)));
}

Functional dWBm2_dn3(double radius) {
  Functional R(radius, "R");
  Functional vtt = VectorThirdTerm(radius);
  Functional n3_ = n3(radius);
  Functional omn3 = 1 - n3_;
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  Functional nV22 = sqr(n2x) + sqr(n2y) + sqr(n2z);
  Functional n22mnV22 = sqr(n2(radius)) - nV22;
  return
    // Psi1
    n0(radius)/omn3
    +
    // Psi2
    n22mnV22/(4*M_PI*R*omn3)
    *((2.0/9.0)*n3_*phi2(n3_) + (1.0/9.0)*sqr_n3_times_phi2_dn3(n3_)) +
    (Functional(1) + (1/9.0)*sqr(n3_)*phi2(n3_)) * (n22mnV22) / (4*M_PI*R*sqr(omn3)) 
    +
    // Psi3
    (Functional(1)/(24*M_PI*sqr(omn3)))*vtt*
    ((-4.0/9.0)*phi3(n3_) - (4.0/9.0)*n3_times_phi3_dn3(n3_)) +
    2*(1 - (4.0/9.0)*n3_*phi3(n3_)) * vtt / (24*M_PI*omn3*omn3*omn3);
}


Functional dAdR_A_over_n(double radius) {
  Functional R(radius, "R");
  Functional two_over_R = Functional(2)/R;
  Functional one_over_4piR = Functional(1)/(4*M_PI*R);
  Functional one_over_4piRsqr = Functional(1)/(4*M_PI*sqr(R));

  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);

  return ShellConvolve(radius)(dWB_dn3(radius))
    + ShellConvolve(radius)(- 2*one_over_4piRsqr/R*dWB_dn0(radius)
                            - one_over_4piRsqr*dWB_dn1(radius))
    + ShellPrimeConvolve(radius)(dWB_dn2(radius)
                                 + one_over_4piR*dWB_dn1(radius)
                                 + one_over_4piRsqr*dWB_dn0(radius))
    + xShellConvolve(radius)(-one_over_4piRsqr*n2x*dWB_dn1v_over_n2v(radius))
    + yShellConvolve(radius)(-one_over_4piRsqr*n2y*dWB_dn1v_over_n2v(radius))
    + zShellConvolve(radius)(-one_over_4piRsqr*n2z*dWB_dn1v_over_n2v(radius))
    + xShellPrimeConvolve(radius)(one_over_4piR*n2x*dWB_dn1v_over_n2v(radius) +
                                  n2x*dWB_dn2v_over_n2v(radius))
    + yShellPrimeConvolve(radius)(one_over_4piR*n2y*dWB_dn1v_over_n2v(radius) +
                                  n2y*dWB_dn2v_over_n2v(radius))
    + zShellPrimeConvolve(radius)(one_over_4piR*n2z*dWB_dn1v_over_n2v(radius) +
                                  n2z*dWB_dn2v_over_n2v(radius));
}

Functional dAdR_S(double radius) {
  Functional R(radius, "R");
  Functional two_over_R = Functional(2)/R;
  Functional one_over_4piR = Functional(1)/(4*M_PI*R);
  Functional one_over_4piRsqr = Functional(1)/(4*M_PI*sqr(R));

  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);

  Functional n2px = xShellPrimeConvolve(radius);
  Functional n2py = yShellPrimeConvolve(radius);
  Functional n2pz = zShellPrimeConvolve(radius);
  Functional n2p = ShellPrimeConvolve(radius);

  return n2(radius)*(dWB_dn3(radius)
                     - dWB_dn1(radius)/(4*M_PI*sqr(R))
                     - dWB_dn0(radius)/(2*M_PI*Pow(3)(R)))
    + n2p*(dWB_dn2(radius)
           + dWB_dn1(radius)/(4*M_PI*R)
           + dWB_dn0(radius)/(4*M_PI*sqr(R)))
    + n2x*(dWB_dn2v_over_n2v(radius)*n2px
           + dWB_dn1v_over_n2v(radius)*n2px/(4*M_PI*R)
           - dWB_dn1v_over_n2v(radius)*n2x/(4*M_PI*sqr(R)))
    + n2y*(dWB_dn2v_over_n2v(radius)*n2py
           + dWB_dn1v_over_n2v(radius)*n2py/(4*M_PI*R)
           - dWB_dn1v_over_n2v(radius)*n2y/(4*M_PI*sqr(R)))
    + n2z*(dWB_dn2v_over_n2v(radius)*n2pz
           + dWB_dn1v_over_n2v(radius)*n2pz/(4*M_PI*R)
           - dWB_dn1v_over_n2v(radius)*n2z/(4*M_PI*sqr(R)));
}


Functional Correlation_A(double radius) {
  Functional R(radius, "R");
  return (Functional(1)/(4*(4*M_PI*sqr(R))))*dAdR_A_over_n(radius)/n_A(radius);
}

Functional Correlation_S(double radius) {
  Functional R(radius, "R");
  return (Functional(1)/(4*(4*M_PI*sqr(R))))*dAdR_S(radius)/sqr(n0(radius));
}

 
Functional dAdR_A_over_n_WBm2(double radius) {
  Functional R(radius, "R");
  Functional two_over_R = Functional(2)/R;
  Functional one_over_4piR = Functional(1)/(4*M_PI*R);
  Functional one_over_4piRsqr = Functional(1)/(4*M_PI*sqr(R));

  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);

  return ShellConvolve(radius)(dWBm2_dn3(radius))
    + ShellConvolve(radius)(- 2*one_over_4piRsqr/R*dWBm2_dn0(radius)
                            - one_over_4piRsqr*dWBm2_dn1(radius))
    + ShellPrimeConvolve(radius)(dWBm2_dn2(radius)
                                 + one_over_4piR*dWBm2_dn1(radius)
                                 + one_over_4piRsqr*dWBm2_dn0(radius))
    + xShellConvolve(radius)(-one_over_4piRsqr*n2x*dWBm2_dn1v_over_n2v(radius))
    + yShellConvolve(radius)(-one_over_4piRsqr*n2y*dWBm2_dn1v_over_n2v(radius))
    + zShellConvolve(radius)(-one_over_4piRsqr*n2z*dWBm2_dn1v_over_n2v(radius))
    + xShellPrimeConvolve(radius)(one_over_4piR*n2x*dWBm2_dn1v_over_n2v(radius) +
                                  n2x*dWBm2_dn2v_over_n2v(radius))
    + yShellPrimeConvolve(radius)(one_over_4piR*n2y*dWBm2_dn1v_over_n2v(radius) +
                                  n2y*dWBm2_dn2v_over_n2v(radius))
    + zShellPrimeConvolve(radius)(one_over_4piR*n2z*dWBm2_dn1v_over_n2v(radius) +
                                  n2z*dWBm2_dn2v_over_n2v(radius));
}

Functional dAdR_S_WBm2(double radius) {
  Functional R(radius, "R");
  Functional two_over_R = Functional(2)/R;
  Functional one_over_4piR = Functional(1)/(4*M_PI*R);
  Functional one_over_4piRsqr = Functional(1)/(4*M_PI*sqr(R));

  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);

  Functional n2px = xShellPrimeConvolve(radius);
  Functional n2py = yShellPrimeConvolve(radius);
  Functional n2pz = zShellPrimeConvolve(radius);
  Functional n2p = ShellPrimeConvolve(radius);

  return n2(radius)*(dWBm2_dn3(radius)
                     - dWBm2_dn1(radius)/(4*M_PI*sqr(R))
                     - dWBm2_dn0(radius)/(2*M_PI*Pow(3)(R)))
    + n2p*(dWBm2_dn2(radius)
           + dWBm2_dn1(radius)/(4*M_PI*R)
           + dWBm2_dn0(radius)/(4*M_PI*sqr(R)))
    + n2x*(dWBm2_dn2v_over_n2v(radius)*n2px
           + dWBm2_dn1v_over_n2v(radius)*n2px/(4*M_PI*R)
           - dWBm2_dn1v_over_n2v(radius)*n2x/(4*M_PI*sqr(R)))
    + n2y*(dWBm2_dn2v_over_n2v(radius)*n2py
           + dWBm2_dn1v_over_n2v(radius)*n2py/(4*M_PI*R)
           - dWBm2_dn1v_over_n2v(radius)*n2y/(4*M_PI*sqr(R)))
    + n2z*(dWBm2_dn2v_over_n2v(radius)*n2pz
           + dWBm2_dn1v_over_n2v(radius)*n2pz/(4*M_PI*R)
           - dWBm2_dn1v_over_n2v(radius)*n2z/(4*M_PI*sqr(R)));
}

Functional Correlation_A_WBm2(double radius) {
  Functional R(radius, "R");
  return (Functional(1)/(4*(4*M_PI*sqr(R))))*dAdR_A_over_n_WBm2(radius)/n_A(radius);
}

Functional Correlation_S_WBm2(double radius) {
  Functional R(radius, "R");
  return (Functional(1)/(4*(4*M_PI*sqr(R))))*dAdR_S_WBm2(radius)/sqr(n0(radius));
}

Functional GrossCorrelation(double radius) {
  Functional R(radius, "R");
  Functional n2prime = ShellConvolve(2*radius, Expression("2*R"));
  Functional n0prime = n2prime/(4*M_PI*sqr(2*R));
  Functional n3prime = StepConvolve(2*radius, Expression("2*R"));

  Functional eta = (1.0/8)*n3prime;
  Functional ghs = (1 - 0.5*eta)/Pow(3)(1 - eta);

  return ghs;
}

Functional YuWuCorrelation_S(double radius) {
  Functional n2 = ShellConvolve(radius);
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  Functional zeta = 1 - (sqr(n2x) + sqr(n2y) + sqr(n2z))/sqr(n2);
  return YuWuCorrelation(radius)*sqr(zeta); // Here we multiply by zeta^2 to get our symmetrical answer.
}

Functional YuWuCorrelation(double radius) {
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

  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  Functional zeta = 1 - (sqr(n2x) + sqr(n2y) + sqr(n2z))/sqr(n2);

  // This gHS (called simply gHS) is the gHS that is used in Fu and
  // Wu's 2005 paper, in equation 13.  It seems ideal, since it
  // includes spatial dependence and is published and tested in
  // various ways.

  // zeta2 is defined right after equation 13 in Fu and Wu 2005.  But
  // it has a typo in it, and should really be the following, which is
  // the packing fraction (and is dimensionless).  Note that this
  // matches the Yu and Wu 2002 paper which is cited by Fu and Wu
  // 2005.
  Functional zeta2 = (R/Functional(3))*n2;
  //Functional invdiff = Functional(1)/(1-zeta3);

  // A careful reading of the Yu and Wu paper indicates that n3 is
  // used rather than n2...
  Functional invdiff = Functional(1)/(1-n3);
  // This is equation 13 in Fu and Wu 2005:
  //return invdiff + 1.5*n3*zeta*sqr(invdiff) + 0.5*sqr(n3)*zeta*Pow(3)(invdiff);

  // This is equation 13 in Fu and Wu 2005, but written to be slightly
  // more efficient:
  Functional ghs = invdiff*(Functional(1) +
                            0.5*(invdiff*zeta2)*zeta*(Functional(3) + invdiff*zeta2));

  return ghs;
}
