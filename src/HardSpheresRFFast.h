// -*- mode: C++; -*-

#pragma once

#include "Functional.h"

class HardSpheresRFFast_type : public FunctionalInterface {
public:
  HardSpheresRFFast_type(double kT_arg, double R_arg) : kT(kT_arg), R(R_arg) { have_analytic_grad = false; }
  double transform(double) const {
    return 0;
  }
  double derive(double) const {
    return 0;
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &x) const {
    assert(&gd); // to avoid an unused parameter error
    assert(&x); // to avoid an unused parameter error
    return 2.7*((-79577.47154594767*ifft(gd, shell(gd, R).cwise()*fft(gd, x))).cwise()*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))).cwise().log()) + 2.7*((ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise().square() - xShellConvolve(R)(gd, x).cwise().square() + (-yShellConvolve(R)(gd, x).cwise().square() - zShellConvolve(R)(gd, x).cwise().square())).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))) + 2.7*((ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise()*(ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise().square() + -3*(xShellConvolve(R)(gd, x).cwise().square() + yShellConvolve(R)(gd, x).cwise().square() + zShellConvolve(R)(gd, x).cwise().square()))).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))).cwise().square()));
  }
  Functional grad(const Functional &, bool) const {
    assert(false);
    return 0;
  }

  void grad(const GridDescription &gd, const VectorXd &x, const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    assert(&gd); // to avoid an unused parameter error
    assert(&x); // to avoid an unused parameter error
    *outgrad += ifft(gd, step(gd, R).cwise()*fft(gd, -(((-79577.47154594767*ifft(gd, shell(gd, R).cwise()*fft(gd, x))).cwise()*(2.7*ingrad)).cwise()/(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x)))))) + (ifft(gd, shell(gd, R).cwise()*fft(gd, -79577.47154594767*((VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))).cwise().log().cwise()*(2.7*ingrad)))) + (2*ifft(gd, shell(gd, R).cwise()*fft(gd, ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise()*((2.7*ingrad).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))))) + (2*xShellConvolve(R)(gd, -(xShellConvolve(R)(gd, x).cwise()*-((2.7*ingrad).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))))) + (2*yShellConvolve(R)(gd, -(yShellConvolve(R)(gd, x).cwise()*-((2.7*ingrad).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))))) + (2*zShellConvolve(R)(gd, -(zShellConvolve(R)(gd, x).cwise()*-((2.7*ingrad).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))))) + (2*ifft(gd, shell(gd, R).cwise()*fft(gd, ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise()*(ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise()*((2.7*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))).cwise().square()))))) - -ifft(gd, step(gd, R).cwise()*fft(gd, -0.01256637061435917*(((ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise().square() - xShellConvolve(R)(gd, x).cwise().square() + (-yShellConvolve(R)(gd, x).cwise().square() - zShellConvolve(R)(gd, x).cwise().square())).cwise()*(2.7*ingrad)).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x)))).cwise().square())))))))));
    if (outpgrad) *outpgrad += ifft(gd, step(gd, R).cwise()*fft(gd, -(((-79577.47154594767*ifft(gd, shell(gd, R).cwise()*fft(gd, x))).cwise()*(2.7*ingrad)).cwise()/(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x)))))) + (ifft(gd, shell(gd, R).cwise()*fft(gd, -79577.47154594767*((VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))).cwise().log().cwise()*(2.7*ingrad)))) + (2*ifft(gd, shell(gd, R).cwise()*fft(gd, ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise()*((2.7*ingrad).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))))) + (2*xShellConvolve(R)(gd, -(xShellConvolve(R)(gd, x).cwise()*-((2.7*ingrad).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))))) + (2*yShellConvolve(R)(gd, -(yShellConvolve(R)(gd, x).cwise()*-((2.7*ingrad).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))))) + (2*zShellConvolve(R)(gd, -(zShellConvolve(R)(gd, x).cwise()*-((2.7*ingrad).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))))))) + (2*ifft(gd, shell(gd, R).cwise()*fft(gd, ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise()*(ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise()*((2.7*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x))).cwise().square()))))) - -ifft(gd, step(gd, R).cwise()*fft(gd, -0.01256637061435917*(((ifft(gd, shell(gd, R).cwise()*fft(gd, x)).cwise().square() - xShellConvolve(R)(gd, x).cwise().square() + (-yShellConvolve(R)(gd, x).cwise().square() - zShellConvolve(R)(gd, x).cwise().square())).cwise()*(2.7*ingrad)).cwise()/(0.01256637061435917*(VectorXd::Ones(gd.NxNyNz) - ifft(gd, step(gd, R).cwise()*fft(gd, x)))).cwise().square())))))))));
  }

  Expression printme(const Expression &) const {
    return Expression("Can't print optimized Functionals");
  }
private:
  double kT;
  double R;
};

inline Functional HardSpheresRFFast(double kT, double R) {
  return Functional(new HardSpheresRFFast_type(kT, R), "HardSpheresRFFast");
}
