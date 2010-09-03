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
#include "Grid.h"

class IdentityType : public FieldFunctionalInterface {
public:
  IdentityType() {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    return data;
  }
  double transform(double n) const {
    return n;
  }

  void grad(const GridDescription &, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += ingrad;
    if (outpgrad) *outpgrad += ingrad; // FIXME: propogate preconditioning
  }
};

FieldFunctional Identity() { return FieldFunctional(new IdentityType()); }

class ChainRuleType : public FieldFunctionalInterface {
public:
  ChainRuleType(const FieldFunctional &fa, const FieldFunctional &fb) : f1(fa), f2(fb) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f1(gd, f2(gd, data));
  }
  double transform(double n) const {
    return f1(f2(n));
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid outgrad1(gd);
    outgrad1.setZero();
    f1.grad(gd, f2(gd, data), ingrad, &outgrad1, 0);
    f2.grad(gd, data, outgrad1, outgrad, outpgrad);
  }
private:
  FieldFunctional f1, f2;
};

FieldFunctional FieldFunctional::operator()(const FieldFunctional &f) const {
  return FieldFunctional(new ChainRuleType(*this, f));
}

class QuotientRuleType : public FieldFunctionalInterface {
public:
  QuotientRuleType(const FieldFunctional &fa, const FieldFunctional &fb) : f1(fa), f2(fb) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f1(gd, data).cwise()/f2(gd, data);
  }
  double transform(double n) const {
    return f1(n)/f2(n);
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    VectorXd out2 = f2(gd, data);
    f1.grad(gd, data, ingrad.cwise()/out2, outgrad, outpgrad);
    f2.grad(gd, data, (ingrad.cwise()*f1(gd, data)).cwise()/((-out2).cwise()*out2), outgrad, outpgrad);
  }
private:
  FieldFunctional f1, f2;
};

FieldFunctional FieldFunctional::operator/(const FieldFunctional &f) const {
  return FieldFunctional(new QuotientRuleType(*this, f));
}

class ProductRuleType : public FieldFunctionalInterface {
public:
  ProductRuleType(const FieldFunctional &fa, const FieldFunctional &fb) : f1(fa), f2(fb) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f1(gd, data).cwise()*f2(gd, data);
  }
  double transform(double n) const {
    return f1(n)*f2(n);
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    f1.grad(gd, data, ingrad.cwise()*f2(gd, data), outgrad, outpgrad);
    f2.grad(gd, data, ingrad.cwise()*f1(gd, data), outgrad, outpgrad);
  }
private:
  FieldFunctional f1, f2;
};

FieldFunctional FieldFunctional::operator*(const FieldFunctional &f) const {
  return FieldFunctional(new ProductRuleType(*this, f));
}

class ScalarRuleType : public FieldFunctionalInterface {
public:
  ScalarRuleType(const FieldFunctional &fa, double xx) : f(fa), x(xx) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return x*f(gd, data);
  }
  double transform(double n) const {
    return x*f(n);
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    f.grad(gd, data, x*ingrad, outgrad, outpgrad);
  }
private:
  FieldFunctional f;
  double x;
};

FieldFunctional FieldFunctional::operator*(double x) const {
  return FieldFunctional(new ScalarRuleType(*this, x));
}


class ScalarMinusType : public FieldFunctionalInterface {
public:
  ScalarMinusType(const FieldFunctional &fa, double xx) : f(fa), x(xx) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return x*VectorXd::Ones(gd.NxNyNz) - f(gd, data);
  }
  double transform(double n) const {
    return x - f(n);
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    f.grad(gd, data, -ingrad, outgrad, outpgrad);
  }
private:
  FieldFunctional f;
  double x;
};

FieldFunctional operator-(double x, const FieldFunctional &f) {
  return FieldFunctional(new ScalarMinusType(f, x));
}


class LogType : public FieldFunctionalInterface {
public:
  LogType(const FieldFunctional &fa) : f(fa) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f(gd, data).cwise().log();
  }
  double transform(double n) const {
    return log(f(n));
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid ingrad1(gd, ingrad);
    ingrad1.cwise() /= f(gd, data);
    f.grad(gd, data, ingrad1, outgrad, outpgrad);
  }
private:
  FieldFunctional f;
};

FieldFunctional log(const FieldFunctional &f) {
  return FieldFunctional(new LogType(f));
}

FieldFunctional FieldFunctional::operator-(const FieldFunctional &f) const {
  return *this + -1*f;
}

class SquareRuleType : public FieldFunctionalInterface {
public:
  SquareRuleType(const FieldFunctional &fa) : f(fa) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f(gd, data).cwise().square();
  }
  double transform(double n) const {
    double x = f(n);
    return x*x;
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    f.grad(gd, data, ingrad.cwise()*(2*f(gd, data)), outgrad, outpgrad);
  }
private:
  FieldFunctional f;
};

FieldFunctional sqr(const FieldFunctional &f) {
  return FieldFunctional(new SquareRuleType(f));
}
