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
#include "new/HomogeneousSFMTFluidFast.h"
#include "equation-of-state.h"
#include "utilities.h"
#include "handymath.h"

// The following tells fac how to run the executable to generate a
// homogeneous.dat file.
/*--

self.rule(exe, [exe], [os.path.dirname(exe)+'/homogeneous.dat'])

--*/

int main(int, char **) {
  FILE *out = fopen("papers/fuzzy-fmt/figs/homogeneous.dat", "w");
  const double Tmax = 10.0, dT = 0.01, Tmin = dT;
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
    for (double T = Tmin; T<= Tmax + dT/2; T += dT) {
      const double temp = T;

      HomogeneousSFMTFluid hf;
      hf.sigma() = 1; // for this computation we will use sigma as our
      // unit of distance to keep things simple.
      hf.epsilon() = 1;
      hf.kT() = temp;
      hf.n() = n_reduced;
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
