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

FieldFunctional HardSpheresRF(double radius, double temperature) {
  FieldFunctional n0 = Identity();
  FieldFunctional n3 = StepConvolve(radius);
  FieldFunctional one_minus_n3 = 1 - StepConvolve(radius);
  FieldFunctional n2 = ShellConvolve(radius);
  FieldFunctional n2x = xShellConvolve(radius);
  FieldFunctional n2y = yShellConvolve(radius);
  FieldFunctional n2z = zShellConvolve(radius);
  FieldFunctional nTxx = xxShellConvolve(radius);
  FieldFunctional nTyy = yyShellConvolve(radius);
  FieldFunctional nTzz = zzShellConvolve(radius);
  FieldFunctional nTxy = xyShellConvolve(radius);
  FieldFunctional nTyz = yzShellConvolve(radius);
  FieldFunctional nTzx = zxShellConvolve(radius*temperature);
  FieldFunctional phi1 = -1*n0*log(one_minus_n3);
  phi1.set_name("phi1");
  const double four_pi_r2 = (4*M_PI)*(radius*radius);
  // n1 is n2/(four_pi_r2)
  FieldFunctional phi2 = (n2*n2 - n2x*n2x - n2y*n2y - n2z*n2z)/(four_pi_r2*one_minus_n3);
  phi2.set_name("phi2");
  FieldFunctional phi3 = n2*(n2*n2 - 3*(n2x*n2x + n2y*n2y + n2z*n2z))/(24*M_PI*one_minus_n3*one_minus_n3);
  phi3.set_name("phi3");
  FieldFunctional total = temperature*(phi1 + phi2 + phi3);
  total.set_name("hard sphere excess");
  return total;
}

FieldFunctional HardSpheresWB(double radius, double temperature) {
  FieldFunctional n0 = Identity();
  FieldFunctional n3 = StepConvolve(radius);
  FieldFunctional one_minus_n3 = 1 - n3;
  FieldFunctional n2 = ShellConvolve(radius);
  FieldFunctional n2x = xShellConvolve(radius);
  FieldFunctional n2y = yShellConvolve(radius);
  FieldFunctional n2z = zShellConvolve(radius);
  FieldFunctional nTxx = xxShellConvolve(radius);
  FieldFunctional nTyy = yyShellConvolve(radius);
  FieldFunctional nTzz = zzShellConvolve(radius);
  FieldFunctional nTxy = xyShellConvolve(radius);
  FieldFunctional nTyz = yzShellConvolve(radius);
  FieldFunctional nTzx = zxShellConvolve(radius*temperature);
  FieldFunctional nTxz = nTzx, nTyx = nTxy, nTzy = nTyz;
  /*
  FieldFunctional trace_nT3 =
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
  FieldFunctional trace_nT3 =
    6*nTxy*nTyz*nTzx +
    nTxx*(  sqr(nTxx) + 3*sqr(nTxy) + 3*sqr(nTzx)) +
    nTyy*(3*sqr(nTxy) +   sqr(nTyy) + 3*sqr(nTyz)) +
    nTzz*(3*sqr(nTzx) + 3*sqr(nTzy) +   sqr(nTzz));
  // * /
  FieldFunctional phi1 = -1*n0*log(one_minus_n3);
  phi1.set_name("phi1");
  const double four_pi_r2 = (4*M_PI)*(radius*radius);
  // n1 is n2/(four_pi_r2)
  FieldFunctional phi2 = (sqr(n2) - sqr(n2x) - sqr(n2y) - sqr(n2z))/(four_pi_r2*one_minus_n3);
  phi2.set_name("phi2");
  FieldFunctional phi3 = (n3 + sqr(one_minus_n3)*log(one_minus_n3))
    /(36*M_PI*sqr(n3)*sqr(one_minus_n3))
    *(n2*sqr(n2) - 3*n2*(sqr(n2x) + sqr(n2y) + sqr(n2z))
      +
      9*(sqr(n2x)*nTxx + sqr(n2y)*nTyy + sqr(n2z)*nTzz
           + 2*(n2x*n2y*nTxy + n2y*n2z*nTyz + n2z*n2x*nTzx)
           - 0.5*trace_nT3));
  phi3.set_name("phi3");
  //FieldFunctional total = temperature*(phi1 + phi2 + phi3);
  //total.set_name("hard sphere excess");
  FieldFunctional total = (temperature*phi1).set_name("phi1") +
    (temperature*phi2).set_name("phi2") + (temperature*phi3).set_name("phi3");
  return total;
}

FieldFunctional HardSpheres(double radius, double temperature) {
  return HardSpheresWB(radius, temperature);
}
