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
#include <math.h>

class GaussianPolynomialType : public FunctionalInterface {
public:
  GaussianPolynomialType(double amplitude, double width, int pow)
    : A(amplitude), sig(width), power(pow) {
    ksig = 1.0/sig; // FIXME: get width right in k space!
  }

  double energy(const GridDescription &gd, const VectorXd &data) const {
    Grid nbar(broaden(gd, data));
    
    double e = 0;
    for (int i=0; i < gd.NxNyNz; i++) {
      const double nbari = nbar[i];
      if (nbari != nbari) {
        printf("nbar[%d] == %g\n", i, nbari);
        assert(nbari == nbari); // check for NaN
      }
      double ehere = A;
      for (int j=0;j<power;j++) ehere *= nbari;
      e += ehere;
    }
    return e*gd.dvolume;
  }
  double energy(double n) const {
    double e = A;
    for (int j=0;j<power;j++) e *= n;
    //printf("\tGaussian energy: %g (n = %g)\n", e, n);
    return e;
  }

  Grid broaden(const GridDescription &gd, const VectorXd &x) const {
    Grid xbar(gd);
    xbar = x;
    ReciprocalGrid xg(xbar.fft());
    xg = xg.cwise()*(-(0.5/ksig/ksig)*xg.g2()).cwise().exp();
    xbar = xg.ifft();
    return xbar;
  }

  void grad(const GridDescription &gd, const VectorXd &n, VectorXd *g_ptr, VectorXd *pg_ptr = 0) const {
    Grid nbar(broaden(gd, n));
    
    Grid dFdnbar(gd);
    for (int i=0; i < gd.NxNyNz; i++) {
      const double nbari = nbar[i];
      dFdnbar[i] = A*power*gd.dvolume;
      for (int j=1; j<power;j++) dFdnbar[i] *= nbari;
    }
    Grid dFdn(broaden(gd, dFdnbar));
    *g_ptr += dFdn;
    if (pg_ptr) *pg_ptr += dFdn;
  }

  void  print_summary(const char *prefix) const {
    printf("%sGaussianPolynomial energy = %g\n", prefix, last_energy);
  }
private:
  double A, sig, ksig;
  int power;
};

Functional GaussianPolynomial(double amplitude, double width, int power) {
  return Functional(new GaussianPolynomialType(amplitude, width, power));
}
