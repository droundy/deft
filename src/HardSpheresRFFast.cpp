#include "Functionals.h"

class HardSpheresRFFast_type : public FunctionalInterface {
public:
  HardSpheresRFFast_type(double R_arg, double kT_arg) : R(R_arg), kT(kT_arg) { have_analytic_grad = false; }
  double transform(double) const {
    return 0;
  }
  double derive(double) const {
    return 0;
  }

  VectorXd transform(const GridDescription &gd, const VectorXd &x) const {
    assert(&gd); // to avoid an unused parameter error
    assert(&x); // to avoid an unused parameter error
    VectorXcd var_0(fft(gd, x));
    VectorXd shell_1(ifft(gd, shell(gd, R).cwise()*var_0));
    VectorXd step_2(ifft(gd, step(gd, R).cwise()*var_0));
    return kT*((-1/(12.56637061435917*(R*R))*shell_1).cwise()*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().log()) + kT*((shell_1.cwise().square() - xShellConvolve(R)(gd, x).cwise().square() + (-yShellConvolve(R)(gd, x).cwise().square() - zShellConvolve(R)(gd, x).cwise().square())).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2))) + kT*((shell_1.cwise()*(shell_1.cwise().square() + -3*(xShellConvolve(R)(gd, x).cwise().square() + yShellConvolve(R)(gd, x).cwise().square() + zShellConvolve(R)(gd, x).cwise().square()))).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square()));
  }
  Functional grad(const Functional &, bool) const {
    assert(false);
    return 0;
  }

  void grad(const GridDescription &gd, const VectorXd &x, const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    assert(&gd); // to avoid an unused parameter error
    assert(&x); // to avoid an unused parameter error
    VectorXcd var_0(fft(gd, x));
    VectorXd shell_1(ifft(gd, shell(gd, R).cwise()*var_0));
    VectorXd step_2(ifft(gd, step(gd, R).cwise()*var_0));
    VectorXcd var_3(fft(gd, -(((-1/(12.56637061435917*(R*R))*shell_1).cwise()*(kT*ingrad)).cwise()/(VectorXd::Ones(gd.NxNyNz) - step_2))));
    VectorXd var_4(ifft(gd, step(gd, R).cwise()*var_3));
    VectorXcd var_5(fft(gd, -1/(12.56637061435917*(R*R))*((VectorXd::Ones(gd.NxNyNz) - step_2).cwise().log().cwise()*(kT*ingrad))));
    VectorXd var_6(ifft(gd, shell(gd, R).cwise()*var_5));
    VectorXcd var_7(fft(gd, shell_1.cwise()*((kT*ingrad).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2)))));
    VectorXd var_8(ifft(gd, shell(gd, R).cwise()*var_7));
    VectorXcd var_9(fft(gd, shell_1.cwise()*(shell_1.cwise()*((kT*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square())))));
    VectorXd var_10(ifft(gd, shell(gd, R).cwise()*var_9));
    VectorXcd var_11(fft(gd, (shell_1.cwise().square() + -3*(xShellConvolve(R)(gd, x).cwise().square() + yShellConvolve(R)(gd, x).cwise().square() + zShellConvolve(R)(gd, x).cwise().square())).cwise()*((kT*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square()))));
    VectorXd var_12(ifft(gd, shell(gd, R).cwise()*var_11));
    VectorXcd var_13(fft(gd, -((VectorXd::Ones(gd.NxNyNz) - step_2).cwise()*(75.39822368615503*(((shell_1.cwise()*(shell_1.cwise().square() + -3*(xShellConvolve(R)(gd, x).cwise().square() + yShellConvolve(R)(gd, x).cwise().square() + zShellConvolve(R)(gd, x).cwise().square()))).cwise()*(kT*ingrad)).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square()).cwise().square())))));
    VectorXd var_14(ifft(gd, step(gd, R).cwise()*var_13));
    VectorXcd var_15(fft(gd, -(12.56637061435917*R*(((shell_1.cwise().square() - xShellConvolve(R)(gd, x).cwise().square() + (-yShellConvolve(R)(gd, x).cwise().square() - zShellConvolve(R)(gd, x).cwise().square())).cwise()*(kT*ingrad)).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2)).cwise().square()))));
    VectorXd var_16(ifft(gd, step(gd, R).cwise()*var_15));
    *outgrad += var_4 + var_6 + (2*var_8 + (2*xShellConvolve(R)(gd, -(xShellConvolve(R)(gd, x).cwise()*-((kT*ingrad).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2))))) + (2*yShellConvolve(R)(gd, -(yShellConvolve(R)(gd, x).cwise()*-((kT*ingrad).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2))))) + (2*zShellConvolve(R)(gd, -(zShellConvolve(R)(gd, x).cwise()*-((kT*ingrad).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2))))) + (2*var_10 + (2*xShellConvolve(R)(gd, -(xShellConvolve(R)(gd, x).cwise()*(3*-(shell_1.cwise()*((kT*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square())))))) + (2*yShellConvolve(R)(gd, -(yShellConvolve(R)(gd, x).cwise()*(3*-(shell_1.cwise()*((kT*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square())))))) + 2*zShellConvolve(R)(gd, -(zShellConvolve(R)(gd, x).cwise()*(3*-(shell_1.cwise()*((kT*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square())))))) + (var_12 + -2*var_14))) - var_16)))));
    if (outpgrad) *outpgrad += var_4 + var_6 + (2*var_8 + (2*xShellConvolve(R)(gd, -(xShellConvolve(R)(gd, x).cwise()*-((kT*ingrad).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2))))) + (2*yShellConvolve(R)(gd, -(yShellConvolve(R)(gd, x).cwise()*-((kT*ingrad).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2))))) + (2*zShellConvolve(R)(gd, -(zShellConvolve(R)(gd, x).cwise()*-((kT*ingrad).cwise()/(12.56637061435917*R*(VectorXd::Ones(gd.NxNyNz) - step_2))))) + (2*var_10 + (2*xShellConvolve(R)(gd, -(xShellConvolve(R)(gd, x).cwise()*(3*-(shell_1.cwise()*((kT*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square())))))) + (2*yShellConvolve(R)(gd, -(yShellConvolve(R)(gd, x).cwise()*(3*-(shell_1.cwise()*((kT*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square())))))) + 2*zShellConvolve(R)(gd, -(zShellConvolve(R)(gd, x).cwise()*(3*-(shell_1.cwise()*((kT*ingrad).cwise()/(75.39822368615503*(VectorXd::Ones(gd.NxNyNz) - step_2).cwise().square())))))) + (var_12 + -2*var_14))) - var_16)))));
  }

  Expression printme(const Expression &) const {
    return Expression("Can't print optimized Functionals");
  }
private:
  double R;
  double kT;
};

Functional HardSpheresRFFast(double R, double kT) {
  return Functional(new HardSpheresRFFast_type(R, kT), "HardSpheresRFFast");
}
