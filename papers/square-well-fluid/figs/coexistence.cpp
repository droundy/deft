// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2015 The Deft Authors
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
#include <cassert>
#include <sys/stat.h>
#include <sys/types.h>

#include "new/HomogeneousSW_liquidFast.h"
#include "handymath.h"

// The following tells fac how to run the executable to generate an
// output file.
/*--

for ww in [1.3]:
    self.rule('%s %g' % (exe,ww), [exe],
              [os.path.dirname(exe)[:-5]+'/data/coexistence/ww%g.dat' % ww])

--*/

const double dT = 0.01;

const double epsilon = 1.0;
const double radius = 1.0;
const double sigma = 2*radius;

const double naccuracy = 1e-5;

enum direction { MINIMIZE, MAXIMIZE };

static inline bool is_improvement(enum direction d, double eold, double enew) {
  if (isnan(eold)) return true;
  if (d == MINIMIZE) {
    return enew < eold;
  } else {
    return enew > eold;
  }
}

static inline double energy(const HomogeneousSW_liquid &f, double n) {
  f.n() = n;
  return f.energy();
}

// Note: we don't assume that nlo < nmid.  Rather, nlo is the fixed
// border, and nmid is a flexible estimate.  It's confusing, but it
// also saves writing two very separate cases.
double optimize(const HomogeneousSW_liquid &hf, enum direction d, double nlo, double nmid) {
  if (d == MINIMIZE) {
    printf("\n\nI am minimizing between %g and %g\n", nlo, nmid);
  } else {
    printf("\n\nI am maximizing between %g and %g\n", nlo, nmid);
  }
  double nhi, emid, elo, ehi;
  elo = energy(hf, nlo);
  emid = energy(hf, nmid);
  if (is_improvement(d, elo, emid)) {
    // great, now we just need to find a high energy!
    nhi = nmid + (nmid - nlo);
    ehi = energy(hf, nhi);
    while (is_improvement(d, emid, ehi)) {
      printf("n is beyond: %10g %10g %10g\n", nlo, nmid, nhi);
      printf("   energies: %10g %10g %10g\n\n", elo, emid, ehi);
      nlo = nmid;
      elo = emid;
      nmid = nhi;
      emid = ehi;
      nhi = 3*nmid - nlo;
      ehi = energy(hf, nhi);
    }
  } else {
    do {
      ehi = emid;
      nhi = nmid;
      printf("n is between: %10g %10g\n", nlo, nhi);
      printf("    energies: %10g %10g\n", elo, ehi);
      nmid = 0.7*nlo + 0.3*nmid;
      emid = energy(hf, nmid);
    } while (is_improvement(d, emid, elo));
  }
  // at this point our minimum should be bracketed between nlo, nhi
  // and nmid, and we reorder lo and hi to be in the right order.
  if (nlo > nhi) {
    double t;
    t = nlo; nlo = nhi; nhi = t;
    t = elo; elo = ehi; ehi = t;
  }
  printf("n bracketed: %10g %10g %10g\n", nlo, nmid, nhi);
  printf("   energies: %10g %10g %10g\n", elo, emid, ehi);
  do {
    if (nmid < 0.5*(nhi+nlo)) {
      const double nnew = 0.6*nhi + 0.4*nlo;
      const double enew = energy(hf, nnew);
      if (is_improvement(d, emid, enew)) {
        elo = emid;
        nlo = nmid;
        emid = enew;
        nmid = nnew;
      } else {
        ehi = enew;
        nhi = nnew;
      }
    } else {
      const double nnew = 0.6*nlo + 0.4*nhi;
      const double enew = energy(hf, nnew);
      if (is_improvement(d, emid, enew)) {
        ehi = emid;
        nhi = nmid;
        emid = enew;
        nmid = nnew;
      } else {
        elo = enew;
        nlo = nnew;
      }
    }
  } while (fabs(nhi - nlo) > naccuracy);
  printf("Final n: %10g\n", 0.5*(nhi+nlo));
  printf(" energy: %10g\n", energy(hf, 0.5*(nhi+nlo)));

  return 0.5*(nhi+nlo);
}

void identify_liquid_vapor(const HomogeneousSW_liquid &hf,
                           double *middlen, // input and output
                           double *nl, double *nv) {
  // first set the slope to zero at middlen:
  hf.n() = *middlen;
  hf.mu() = 0;
  hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf
  double emid, ev, el;
  const double tol = 1e-8;
  do {
    // Find estimates for the liquid and vapor densities.
    *nv = optimize(hf, MINIMIZE, *middlen, 0);
    ev = energy(hf, *nv);
    printf("  got nv: %10g\n", *nv);
    printf("  energy: %10g\n", ev);
    *nl = optimize(hf, MINIMIZE, *middlen, max(2*(*nl), 3*(*middlen)));
    el = energy(hf, *nl);
    printf("  got nl: %10g\n", *nl);
    printf("  energy: %10g\n", el);
    // Adjust the chemical potential to make the free energies at the
    // two minima equal.
    hf.mu() += (ev - el)/(*nv - *nl);
    // Now we find the new middle point that maximizes the free
    // energy.
    *middlen = optimize(hf, MAXIMIZE, *nl, *middlen);
    emid = energy(hf, *middlen);
    printf("  got nm: %10g\n", *middlen);
    printf("  energy: %10g\n", emid);

    // for (double nn=*nv/2; nn < 1.1*(*nl); nn *= 1.2) {
    //   printf("    %10g %10g\n", nn, energy(hf, nn));
    // }

  } while (fabs(ev - el)/fabs(emid-el) > tol);
}

int main(int argc, char **argv) {
  double lambda;
  if (argc != 2) {
    printf("usage: %s lambda\n", argv[0]);
    return 1;
  }

  sscanf(argv[1], "%lg", &lambda);

  printf("lam %g\n", lambda);

  char *fname = new char[4096];
  mkdir("papers/square-well-fluid/data/coexistence", 0777);
  sprintf(fname, "papers/square-well-fluid/data/coexistence/ww%g.dat", lambda);
  FILE *out = fopen(fname, "w");
  assert(out);

  double eta_middle = 0.3;
  double n_middle = eta_middle/(M_PI*uipow(sigma, 3)/6);
  for (double temp = dT; temp < 5; temp += dT) {
    HomogeneousSW_liquid hf;
    hf.R() = radius;
    hf.epsilon() = epsilon;
    hf.kT() = temp;
    hf.lambda() = lambda;
    double nl, nv;
    identify_liquid_vapor(hf, &n_middle, &nl, &nv);
    printf("\n\n**********************************************\n");
    printf("kT %g nl %g nm %g nv %g\n", temp, nl, n_middle, nv);
    printf("**********************************************\n");
    if (isnan(n_middle) || n_middle < nv || n_middle > nl) {
      printf("Got craziness, I am quitting now.\n");
      exit(0);
    }
    fprintf(out, "%g\t%g\t%g\t%g\n", temp, nl, n_middle, nv);
    fflush(out);

    sprintf(fname, "papers/square-well-fluid/data/coexistence/ww%g-kT%g.dat",
            lambda, temp);
    FILE *ff = fopen(fname, "w");
    assert(ff);
    for (double nn=0.1*nv; nn<1.1*nl; nn *= 1.01) {
      fprintf(ff, "%g\t%g\n", nn, energy(hf, nn));
    }
    fclose(ff);
  }

  fclose(out);
  delete[] fname;
  return 0;
}
