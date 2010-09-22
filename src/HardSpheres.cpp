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
#include <stdio.h>
#include <math.h>

Functional HardSpheresRF(double radius, double temperature) {
  Functional R(radius);
  R.set_name("R");
  const Functional four_pi_r = (4*M_PI)*R;
  const Functional four_pi_r2 = (4*M_PI)*sqr(R);
  Functional n3 = StepConvolve(radius);
  Functional one_minus_n3 = 1 - StepConvolve(radius);
  Functional n2 = ShellConvolve(radius);
  Functional n2x = xShellConvolve(radius);
  Functional n2y = yShellConvolve(radius);
  Functional n2z = zShellConvolve(radius);
  Functional phi1 = (Functional(-1)/four_pi_r2)*n2*log(one_minus_n3);
  Functional kT(temperature);
  kT.set_name("kT");
  phi1.set_name("phi1");
  // n1 is n2/(four_pi_r2)
  Functional phi2 = (sqr(n2) - sqr(n2x) - sqr(n2y) - sqr(n2z))/(four_pi_r*one_minus_n3);
  phi2.set_name("phi2");
  Functional phi3 = n2*(sqr(n2) - 3*(sqr(n2x) + sqr(n2y) + sqr(n2z)))/(24*M_PI*sqr(one_minus_n3));
  phi3.set_name("phi3");
  //Functional total = temperature*(phi1 + phi2 + phi3);
  Functional total = (kT*phi1).set_name("phi1") +
    (kT*phi2).set_name("phi2") + (kT*phi3).set_name("phi3");
  return total;
}

Functional HardSpheresWB(double radius, double temperature) {
  Functional R(radius);
  R.set_name("R");
  const Functional four_pi_r = (4*M_PI)*R;
  const Functional four_pi_r2 = (4*M_PI)*sqr(R);
  Functional n3 = StepConvolve(radius);
  Functional one_minus_n3 = 1 - n3;
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
  Functional kT(temperature);
  kT.set_name("kT");
  Functional trace_nT3 =
    6*nTxy*nTyz*nTzx +
    nTxx*(  sqr(nTxx) + 3*sqr(nTxy) + 3*sqr(nTzx)) +
    nTyy*(3*sqr(nTxy) +   sqr(nTyy) + 3*sqr(nTyz)) +
    nTzz*(3*sqr(nTzx) + 3*sqr(nTzy) +   sqr(nTzz));
  // */
  // n0 is n2/(four_pi_r2)
  Functional phi1 = (Functional(-1)/four_pi_r2)*n2*log(one_minus_n3);
  phi1.set_name("phi1");
  // n1 is n2/(four_pi_r)
  Functional phi2 = (sqr(n2) - sqr(n2x) - sqr(n2y) - sqr(n2z))/(four_pi_r*one_minus_n3);
  phi2.set_name("phi2");
  Functional phi3 = (n3 + sqr(one_minus_n3)*log(one_minus_n3))
    /(36*M_PI*sqr(n3)*sqr(one_minus_n3))
    *(n2*sqr(n2) - 3*n2*(sqr(n2x) + sqr(n2y) + sqr(n2z))
      +
      9*(sqr(n2x)*nTxx + sqr(n2y)*nTyy + sqr(n2z)*nTzz
           + 2*(n2x*n2y*nTxy + n2y*n2z*nTyz + n2z*n2x*nTzx)
           - 0.5*trace_nT3));
  phi3.set_name("phi3");
  //Functional total = kT*(phi1 + phi2 + phi3);
  //total.set_name("hard sphere excess");
  Functional total = (kT*phi1).set_name("phi1") +
    (kT*phi2).set_name("phi2") + (kT*phi3).set_name("phi3");
  return total;
}

Functional HardSpheres(double radius, double temperature) {
  return HardSpheresWB(radius, temperature);
}
