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
    if (outpgrad) {
      Grid outgrad1(gd);
      Grid outpgrad1(gd);
      VectorXd out2 = f2(gd, data);
      outgrad1.setZero();
      outpgrad1.setZero();
      f1.grad(gd, data, ingrad, &outgrad1, &outpgrad1);
      *outgrad += out2.cwise()*outgrad1;
      *outpgrad += out2.cwise()*outpgrad1;

      outgrad1.setZero();
      outpgrad1.setZero();
      out2 = f1(gd, data);
      f2.grad(gd, data, ingrad, &outgrad1, &outpgrad1);
      *outgrad += out2.cwise()*outgrad1;
      *outpgrad += out2.cwise()*outpgrad1;
    } else {
      Grid outgrad1(gd);
      VectorXd out2 = f2(gd, data);
      outgrad1.setZero();
      f1.grad(gd, data, ingrad, &outgrad1, 0);
      *outgrad += out2.cwise()*outgrad1;

      outgrad1.setZero();
      out2 = f1(gd, data);
      f2.grad(gd, data, ingrad, &outgrad1, 0);
      *outgrad += out2.cwise()*outgrad1;
    }
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
    if (outpgrad) {
      Grid outgrad1(gd);
      Grid outpgrad1(gd);
      outgrad1.setZero();
      outpgrad1.setZero();
      f.grad(gd, data, ingrad, &outgrad1, &outpgrad1);
      *outgrad += x*outgrad1;
      *outpgrad += x*outpgrad1;
    } else {
      Grid outgrad1(gd);
      outgrad1.setZero();
      f.grad(gd, data, ingrad, &outgrad1, 0);
      *outgrad += x*outgrad1;
    }
  }
private:
  FieldFunctional f;
  double x;
};

FieldFunctional FieldFunctional::operator*(double x) const {
  return FieldFunctional(new ScalarRuleType(*this, x));
}
