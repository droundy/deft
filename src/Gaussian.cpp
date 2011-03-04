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

class GaussianType : public FunctionalInterface {
public:
  GaussianType(double w) {
    width = w;
    kfac = -0.5*width*width; // FIXME: get width right in k space!
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &, const VectorXd &data) const {
    Grid out(gd, data);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= (kfac*g2(gd)).cwise().exp();
    return recip.ifft();
  }
  double transform(double, double n) const {
    return n;
  }
  double derive(double, double) const {
    return 1;
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Functional grad(const Functional &ingrad, const Functional &, bool) const {
    return Gaussian(width)(ingrad);
  }
  Functional grad_T(const Functional &) const {
    return 0;
  }
  void grad(const GridDescription &gd, const VectorXd &, const VectorXd &,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid out(gd, ingrad);
    ReciprocalGrid recip = out.fft();
    recip.cwise() *= (kfac*g2(gd)).cwise().exp();
    out = recip.ifft();
    *outgrad += out;

    // FIXME: we will want to propogate preexisting preconditioning
    if (outpgrad) *outpgrad += out;
  }
  Expression printme(const Expression &x) const {
    Expression out = funexpr("Gaussian", Expression("width"))(Expression("gd"), x);
    out.unlazy = true;
    return out;
  }
  Expression cwiseprintme(const Expression &x) const {
    return x;
  }
private:
  double width, kfac;
};

Functional Gaussian(double width) {
  return Functional(new GaussianType(width));
}

Functional StepConvolve(double R, Expression r) {
  r.set_type("double");
  return Functional(function_for_convolve<step_op<complex> >, R,
                    r, Expression(4*M_PI/3)*(r*r*r), true);
}

Functional ShellConvolve(double R, Expression r) {
  r.set_type("double");
  return Functional(function_for_convolve<shell_op<complex> >, R,
                    r, Expression(4*M_PI)*(r*r), true);
}

Functional xShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<xshell_op<complex> >, R, r, Expression(0), false);
}
Functional yShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<yshell_op<complex> >, R, r, Expression(0), false);
}
Functional zShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<zshell_op<complex> >, R, r, Expression(0), false);
}

Functional xyShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<xyshell_op<complex> >, R, r, Expression(0), true);
}
Functional yzShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<yzshell_op<complex> >, R, r, Expression(0), true);
}
Functional zxShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<zxshell_op<complex> >, R, r, Expression(0), true);
}

Functional xxShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<xxshell_op<complex> >, R, r, Expression(0), true);
}
Functional yyShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<yyshell_op<complex> >, R, r, Expression(0), true);
}
Functional zzShellConvolve(double R, Expression r) {
  return Functional(function_for_convolve<zzshell_op<complex> >, R, r, Expression(0), true);
}
