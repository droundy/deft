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

class PowType : public FunctionalInterface {
public:
  PowType(int nn) : n(nn) {}

  VectorXd transform(double, const GridDescription &gd, const VectorXd &data) const {
    switch (n) {
    case 0: return VectorXd::Ones(data.rows());
    case 1: return data;
    }
    VectorXd out(data);
    for (int i=0; i<gd.NxNyNz; i++)
      for (int p=1; p < n; p++)
        out[i] *= data[i];
    return out;
  }
  double transform(double, double x) const {
    switch (n) {
    case 0: return 1;
    case 1: return x;
    }
    double v = x;
    for (int p=1; p < n; p++) v *= x;
    return v;
  }
  double derive(double, double x) const {
    if (n < 1) return 0;
    double v = n;
    for (int p=1; p < n; p++) v *= x;
    return v;
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool) const {
    switch (n) {
    case 0: return 0;
    case 1: return ingrad;
    }
    return Pow(n-1)(x)*n*ingrad;
  }
  void grad(double, const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    switch (n) {
    case 0: return; // zero gradient!
    case 1:
      *outgrad += ingrad;
      if (outpgrad) *outpgrad += ingrad;
      return;
    }
    for (int i=0; i<gd.NxNyNz; i++) {
      double foo = n*ingrad[i];
      for (int p=1; p < n; p++) foo *= data[i];
      (*outgrad)[i] += foo;
      if (outpgrad) (*outpgrad)[i] += foo;
    }
  }
  Expression printme(const Expression &x) const {
    switch (n) {
    case 0: return 0;
    case 1: return x;
    case 2: return sqr(x);
    }
    // This is more than a little hokey...
    if (n & 1) return x*Pow(n-1).printme(x);
    return Pow(n/2).printme(sqr(x));
  }
private:
  int n;
};

Functional Pow(int n) {
  assert(n > 0);
  return Functional(new PowType(n));
}


class PowAndHalfType : public FunctionalInterface {
public:
  PowAndHalfType(int nn) : n(nn) {}

  VectorXd transform(double, const GridDescription &, const VectorXd &data) const {
    VectorXd out(data.cwise().sqrt());
    if (n < 0) {
      for (int p=0; p > n; p--) out = out.cwise() / data;
    } else {
      for (int p=0; p < n; p++) out = out.cwise() * data;
    }
    return out;
  }
  double transform(double, double x) const {
    double out = sqrt(x);
    if (n < 0) {
      for (int p=0; p > n; p--) out /= x;
    } else {
      for (int p=0; p < n; p++) out *= x;
    }
    return out;
  }
  double derive(double, double x) const {
    double out = (n+0.5)/sqrt(x);
    if (n < 0) {
      for (int p=0; p > n; p--) out /= x;
    } else {
      for (int p=0; p < n; p++) out *= x;
    }
    return out;
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool) const {
    return PowAndHalf(n-1)(x)*(n+0.5)*ingrad;
  }
  void grad(double, const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    if (n > 0) {
      for (int i=0; i<gd.NxNyNz; i++) {
        double foo = (n+0.5)*ingrad[i]/sqrt(data[i]);
        for (int p=1; p < n; p++) foo *= data[i];
        (*outgrad)[i] += foo;
        if (outpgrad) (*outpgrad)[i] += foo;
      }
    } else {
      for (int i=0; i<gd.NxNyNz; i++) {
        double di = data[i];
        double oodi = 1.0/di; // one over d_i
        double foo = (n+0.5)*ingrad[i]/sqrt(di);
        for (int p=0; p > n; p--) foo *= oodi;
        (*outgrad)[i] += foo;
        if (outpgrad) (*outpgrad)[i] += foo;
      }
    }
  }
  Expression printme(const Expression &x) const {
    if (n == 0) {
      return sqrt(x);
    } else if (n >= 0) {
      return Pow(n).printme(x)*sqrt(x);
    } else {
      return sqrt(x) / Pow(-n).printme(x);
    }
  }
private:
  int n;
};

Functional PowAndHalf(int n) {
  return Functional(new PowAndHalfType(n));
}

Functional sqrt(const Functional &f) {
  return PowAndHalf(0)(f);
}
