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

#include <stdio.h>
#include <time.h>
#include "new/HomogeneousWhiteBearFluidFast.h"
#include "equation-of-state.h"
#include "utilities.h"
#include "handymath.h"

// The following tells fac how to run the executable to generate a
// homogeneous-bh.dat file.
/*--

self.rule(exe, [exe], [os.path.dirname(exe)+'/homogeneous-bh.dat'])

--*/

const double eps = 1.0;
const double sigma = 1.0; /* Let's define sigma == 1 for this one? */
const double R = sigma*pow(2.0,1.0/6.0);
const int N = 10000;

double R_BH(const double kT) {
  double bh_diameter = 0;
  const double dr = R/N;
  const double beta = 1.0/kT;
  for (double r_cur=dr/2; r_cur < R; r_cur += dr) {
    const double s6 = uipow(sigma/r_cur, 6);
    bh_diameter += (1 - exp(-beta*(4*eps*(s6*s6 - s6) + eps)))*dr;
  }
  return bh_diameter/2;
}

int main(int, char **) {
  FILE *out = fopen("papers/fuzzy-fmt/figs/homogeneous-bh.dat", "w");
  const double Tmax = 3.1, dT = 0.01, Tmin = dT;
  fprintf(out, "# n_reduced");
  for (double T = Tmin; T<= Tmax + dT/2; T += dT) {
    fprintf(out, "\tp(kT=%g)/nkT", T);
  }
  fprintf(out, "\n0");
  for (double T = Tmin; T<= Tmax + dT/2; T += dT) {
    fprintf(out, "\t%g", T);
  }
  fprintf(out, "\n");

  const double dn = 0.01, nmax = 2.5;
  for (double n_reduced = dn; n_reduced <= nmax; n_reduced += dn) {
    fprintf(out, "%g", n_reduced);
    printf("n_reduced == %g\n", n_reduced);
    for (double T = Tmin; T<= Tmax + dT/2; T += dT) {
      const double temp = T;

      HomogeneousWhiteBearFluid hf;

      hf.R() = R_BH(temp);
      //printf("rad is %g\n", hf.R());
      hf.kT() = temp;
      hf.n() = n_reduced/uipow(sigma,3);
      //printf("dividing by sigma = %g\n", sigma);
      //printf("eta is %g\n", hf.n()*uipow(hf.R(),3)*M_PI*4/3);
      hf.mu() = 0;
      hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf

      const double pressure = n_reduced*hf.d_by_dn() - hf.energy();

      // return *reduced* pressure!
      fprintf(out, "\t%g", pressure);
    }
    fprintf(out, "\n");
  }
  fclose(out);
  return 0;
}
