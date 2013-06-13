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
#include "handymath.h"
#include "Grid.h"

bool FunctionalInterface::I_have_analytic_grad() const {
  return true;
}

bool FunctionalInterface::I_am_homogeneous() const {
  return false;
}

bool FunctionalInterface::I_am_local() const {
  return true;
}

bool FunctionalInterface::I_am_zero() const {
  return false;
}

bool FunctionalInterface::I_am_one() const {
  return false;
}

bool FunctionalInterface::I_am_constant_wrt_x() const {
  return false;
}

bool FunctionalInterface::I_preserve_homogeneous() const {
  return true;
}

bool FunctionalInterface::I_give_zero_for_zero() const {
  return false;
}

void FunctionalInterface::pgrad(const GridDescription &gd, double kT, const VectorXd &x,
                                const VectorXd &ingrad,
                                VectorXd *outpgrad) const {
  Grid trash(gd);
  grad(gd, kT, x, ingrad, &trash, outpgrad);
}

void FunctionalInterface::print_summary(const char *prefix, double e, std::string name) const {
  if (name != "") printf("%s%25s =", prefix, name.c_str());
  else printf("%s%25s =", prefix, "UNKNOWN");
  print_double("", e);
  printf("\n");
}

double FunctionalInterface::integral(const GridDescription &gd, double kT, const VectorXd &x) const {
  return transform(gd, kT, x).sum()*gd.dvolume;
}

// The following is a "fake" functional, used for dumping code to
// generate the gradient.
class PretendIngradType : public FunctionalInterface {
public:
  PretendIngradType() {}

  VectorXd transform(const GridDescription &, double, const VectorXd &x) const {
    return x;
  }
  double transform(double, double) const {
    return 1;
  }
  double derive(double, double) const {
    return 1; // hokey!
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return Functional(double(0));
  }
  Functional grad_T(const Functional &) const {
    return Functional(0.0);
  }
  void grad(const GridDescription &, double, const VectorXd &, const VectorXd &,
            VectorXd *, VectorXd *) const {
  }
  bool I_am_homogeneous() const {
    return true; // FIXME:  this should be false, but our code would be too slow!
  }
};


int Functional::run_finite_difference_test(const char *testname, double temp, const Grid &x,
                                           const VectorXd *direction) const {
  printf("\nRunning finite difference test on %s:\n", testname);
  const GridDescription gd(x.description());

  int retval = 0;

  VectorXd my_grad(x);
  double Eold = integral(temp, x);
  my_grad.setZero();
  integralgrad(temp, x, &my_grad);
  if (my_grad.norm() == 0.0) {
    printf("FAIL: Gradient is zero...\n");
    return 1;
  }
  VectorXd my_direction(my_grad);
  if (direction) my_direction = *direction;
  my_direction /= my_direction.norm();

  const double lderiv = my_direction.dot(my_grad);
  if (lderiv == 0.0) {
    printf("FAIL: Gradient is zero...\n");
    retval++;
  }
  {
    // compare the gradient computed by integralgrad() with and
    // without a non-null pgrad argument
    VectorXd my_pgrad(x);
    my_grad.setZero();
    my_pgrad.setZero();
    integralgrad(temp, x, &my_grad, &my_pgrad);
    const double lderiv_new = my_direction.dot(my_grad);
    if (fabs(lderiv_new/lderiv - 1) > 1e-12) {
      printf("\n*** WARNING!!! INCONSISTENT GRADIENTS! ***\n");
      printf("Different gradient in integralgrad() with and without pgradout\n");
      printf("FAIL: Fractional error is %g\n\n", lderiv_new/lderiv - 1);
      retval++;
    }

    if (I_have_analytic_grad()) {
      // We have an analytic grad, so let's make sure it matches the
      // other one...
      VectorXd othergrad(grad(dV(), Identity(), false)(temp, x));
      double maxgraderr = (othergrad-my_grad).cwise().abs().maxCoeff();
      double maxgrad = my_grad.cwise().abs().maxCoeff();
      if (maxgraderr/maxgrad > 1e-9) {
        printf("othergrad[0] = %g\n", othergrad[0]);
        printf("my_grad[0] = %g\n", my_grad[0]);
        printf("maxgraderr is %g while maxgrad is %g\n", maxgraderr, maxgrad);
        printf("FAIL: Discrepancy in my gradient is too big: %g\n\n", maxgraderr/maxgrad);
        retval++;
      }
    }
  }

  // We try to choose as small values for epsilon as are consistent with
  // getting meaningful results:
  const double sigfigs_best = 1e-15;
  const double epsilon_max = 1e-15*fabs(Eold)/fabs(lderiv)/sigfigs_best;
  const int min_p = (int) -floor(log10(epsilon_max) + 0.5);
  const int max_p = min_p + 7;
  printf("choosing delta of 1e%d based on directional grad %g and energy %g\n",
         -min_p, lderiv, Eold);
  double min = 1e300, best_ratio_error = 1.0;
  VectorXd grads(max_p + 1 - min_p);
  VectorXd diff_grads(max_p - min_p);
  // Take steps in the grad direction.
  for(int p=min_p; p <= max_p; p++) {
    const double eps_ratio = 10.0;
    const double epsilon = pow(eps_ratio, -p);
    // The following is a little wasteful of memory...
    const double Eplus= integral(temp, gd, x + epsilon*my_direction);
    const double Eminus=integral(temp, gd, x - epsilon*my_direction);

    grads[p-min_p] = (Eplus-Eminus)/(2*epsilon);
    printf("    eps^2 = %25.16f deltaE %.12g\n",
           epsilon*epsilon*pow(eps_ratio, 2.0*min_p), Eplus-Eminus);
    printf("FD   Ratio: %25.16f (grad ~ %g)\n", grads[p-min_p]/lderiv,
           grads[p-min_p]);
    printf("FD sigfigs: %25.16f\n", 1e-15*fabs(Eold/(Eplus-Eminus)));
    fflush(stdout);

    if (fabs(grads[p-min_p]/lderiv - 1.0) < best_ratio_error)
      best_ratio_error = fabs(grads[p-min_p]/lderiv - 1.0);
    if (p > min_p) {
      diff_grads[p-min_p-1]=(grads[p-min_p]-grads[p-min_p-1]);
      double diff = (grads[p-min_p-1]-lderiv)*
        (-1+1./(eps_ratio*eps_ratio))/diff_grads[p-min_p-1];
      if ( min > fabs(diff-1) ) min = fabs(diff-1);
      //dft_log("FD diff: %25.16lf\n\n", diff);
    }
  }

  if (min < 1e-3 && best_ratio_error < 1e-7) {
    printf("Passed on basis of reasonable scaling (%g) and accuracy (%g).\n",
           min, best_ratio_error);
  } else if (best_ratio_error < 1e-10) {
    printf("Passed on basis of a gradient accuracy of (%g) (with scaling %g).\n",
           best_ratio_error, min);
  } else if (min < 1e-5 && best_ratio_error < 0.01) {
    printf("Passed on basis of seriously nice scaling %g (and low accuracy %g).\n",
           min, best_ratio_error);
  } else {
    printf("FAIL: Failed with scaling ratio of %g and gradient accuracy of %g\n",
           min, best_ratio_error);
    retval++;
  }

  return retval;
}

int Functional::run_homogeneous_finite_difference_test(const char *testname, double temp, double n) const {
  printf("\nRunning homogeneous finite difference test on %s:\n", testname);
  fflush(stdout);
  int retval = 0;
  const double Eold = (*this)(temp, n);
  const double lderiv = derive(temp, n);

  // We try to choose as small values for epsilon as are consistent with
  // getting meaningful results:
  const double sigfigs_best = 1e-15;
  const double epsilon_max = 1e-15*fabs(Eold)/fabs(lderiv)/sigfigs_best;
  const int min_p = (int) -floor(log10(epsilon_max) + 0.5);
  const int max_p = min_p + 7;
  if (lderiv == 0) {
    printf("skipping finite difference test with zero derivative...\n");
    return 0;
  }
  printf("choosing delta of 1e%d based on derivative %g and energy %g at x = %g\n",
         -min_p, lderiv, Eold, n);
  double min = 1e300, best_ratio_error = 1.0;
  VectorXd grads(max_p + 1 - min_p);
  VectorXd diff_grads(max_p - min_p);
  // Take steps in the grad direction.
  for(int p=min_p; p <= max_p; p++) {
    const double eps_ratio = 10.0;
    const double epsilon = pow(eps_ratio, -p);
    // The following is a little wasteful of memory...
    const double Eplus= (*this)(temp, n + epsilon);
    const double Eminus=(*this)(temp, n - epsilon);

    grads[p-min_p] = (Eplus-Eminus)/(2*epsilon);
    printf("    eps^2 = %25.16f deltaE %.12g\n",
           epsilon*epsilon*pow(eps_ratio, 2.0*min_p), Eplus-Eminus);
    printf("fd   Ratio: %25.16f (grad ~ %g)\n", grads[p-min_p]/lderiv,
           grads[p-min_p]);
    printf("fd sigfigs: %25.16f\n", 1e-15*fabs(Eold/(Eplus-Eminus)));
    fflush(stdout);

    if (fabs(grads[p-min_p]/lderiv - 1.0) < best_ratio_error)
      best_ratio_error = fabs(grads[p-min_p]/lderiv - 1.0);
    if (p > min_p) {
      diff_grads[p-min_p-1]=(grads[p-min_p]-grads[p-min_p-1]);
      double diff = (grads[p-min_p-1]-lderiv)*
        (-1+1./(eps_ratio*eps_ratio))/diff_grads[p-min_p-1];
      if ( min > fabs(diff-1) ) min = fabs(diff-1);
      //dft_log("FD diff: %25.16lf\n\n", diff);
    }
  }

  if (min < 1e-3 && best_ratio_error < 1e-7) {
    printf("Passed on basis of reasonable scaling (%g) and accuracy (%g).\n",
           min, best_ratio_error);
  } else if (best_ratio_error < 1e-10) {
    printf("Passed on basis of a gradient accuracy of (%g) (with scaling %g).\n",
           best_ratio_error, min);
  } else if (min < 1e-5 && best_ratio_error < 0.01) {
    printf("Passed on basis of seriously nice scaling %g (and low accuracy %g).\n",
           min, best_ratio_error);
  } else {
    printf("FAIL: Failed with scaling ratio of %g and gradient accuracy of %g\n",
           min, best_ratio_error);
    retval++;
  }
  return retval;
}

double Functional::print_iteration(const char *prefix, int iter) const {
  printf("%s==============\n", prefix);
  printf("%sIteration %4d\n", prefix, iter);
  printf("%s==============\n", prefix);

  print_summary(prefix, itsCounter->last_energy, get_name());
  const Functional *nxt = this;
  double etot = 0;
  while (nxt) {
    etot += nxt->itsCounter->last_energy;
    //printf("In print_iteration summarizing %p\n", nxt);
    //nxt->print_summary(prefix, nxt->itsCounter->last_energy);
    nxt = nxt->next();
  }
  return etot;
}

void Functional::print_summary(const char *prefix, double energy, std::string name) const {
  if (get_name() != "") name = get_name();
  if (!next()) {
    itsCounter->ptr->print_summary(prefix, energy, name);
    return;
  }
  const Functional *nxt = this;
  double etot = 0;
  while (nxt) {
    if (nxt->get_name() != "") name = nxt->get_name();
    etot += nxt->itsCounter->last_energy;
    nxt->itsCounter->ptr->print_summary(prefix, nxt->itsCounter->last_energy, name);
    nxt = nxt->next();
  }
  printf("%s%25s =", prefix, "total energy");
  print_double("", etot);
  printf("\n");
  //assert(fabs(etot - energy) < 1e-6);
}

Functional Identity() { return Pow(1); }

class dVType : public FunctionalInterface {
public:
  dVType() {}

  VectorXd transform(const GridDescription &gd, double, const VectorXd &) const {
    return gd.dvolume*VectorXd::Ones(gd.NxNyNz);
  }
  double transform(double, double) const {
    return 1;
  }
  double derive(double, double) const {
    return 1; // hokey!
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return Functional(0.0);
  }
  Functional grad_T(const Functional &) const {
    return Functional(0.0);
  }
  void grad(const GridDescription &, double, const VectorXd &, const VectorXd &,
            VectorXd *, VectorXd *) const {
  }
  bool I_am_homogeneous() const {
    return true;
  }
};

Functional dV() {
  return Functional(new dVType());
}

class Constant : public FunctionalInterface {
public:
  Constant(double x, const char *n) : c(x), name(n) {}

  VectorXd transform(const GridDescription &, double, const VectorXd &data) const {
    return c*VectorXd::Ones(data.rows());
  }
  double transform(double, double) const {
    return c;
  }
  double derive(double, double) const {
    return 0;
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return Functional(0.0);
  }
  Functional grad_T(const Functional &) const {
    return Functional(0.0);
  }
  void grad(const GridDescription &, double, const VectorXd &, const VectorXd &,
            VectorXd *, VectorXd *) const {
  }
  bool I_am_homogeneous() const {
    return true;
  }
  bool I_am_zero() const {
    return c == 0 && !name;
  }
  bool I_am_one() const {
    return c == 1 && !name;
  }
  bool I_am_constant_wrt_x() const {
    return true;
  }
private:
  double c;
  const char *name;
};

Functional::Functional(double x, const char *n) : itsCounter(0) {
  init(new Constant(x, n), "");
}

class ConstantField : public FunctionalInterface {
public:
  ConstantField(const VectorXd &x) : c(x) {}
  bool I_preserve_homogeneous() const {
    return false;
  }
  bool I_am_constant_wrt_x() const {
    return true;
  }

  VectorXd transform(const GridDescription &, double, const VectorXd &) const {
    return c;
  }
  double transform(double, double) const {
    return c.sum()/c.rows();
  }
  double derive(double, double) const {
    return 0;
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Functional grad(const Functional &, const Functional &, bool) const {
    return Functional(double(0));
  }
  Functional grad_T(const Functional &) const {
    return Functional(0.0);
  }
  void grad(const GridDescription &, double, const VectorXd &, const VectorXd &,
            VectorXd *, VectorXd *) const {
  }
private:
  VectorXd c;
};

Functional::Functional(const VectorXd &x) : itsCounter(0) {
  init(new ConstantField(x), "");
}

class ChainRuleType : public FunctionalInterface {
public:
  ChainRuleType(const Functional &fa, const Functional &fb) : f1(fa), f2(fb) {}
  bool I_am_local() const {
    return f1.I_am_local() && f2.I_am_local();
  }
  bool I_preserve_homogeneous() const {
    return f1.I_preserve_homogeneous() && f2.I_preserve_homogeneous();
  }
  bool I_am_constant_wrt_x() const {
    return f1.I_preserve_homogeneous() && f2.I_am_constant_wrt_x();
  }

  double integral(const GridDescription &gd, double kT, const VectorXd &data) const {
    if (!f1.next()) {
      // This is the simple, efficient case!
      return f1.integral(gd, kT, f2(gd, kT, data));
    }
    // This does some extra work to save the energies of each term in
    // the sum, just in case we want to print it!
    VectorXd f2data(f2(gd, kT, data));
    double e = f1.justMeIntegral(gd, kT, f2data);
    f1.set_last_energy(e);
    Functional *nxt = f1.next();
    while (nxt) {
      double enxt = nxt->justMeIntegral(gd, kT, f2data);
      nxt->set_last_energy(enxt);
      e += enxt;
      nxt = nxt->next();
    }
    return e;
  }
  VectorXd transform(const GridDescription &gd, double kT, const VectorXd &data) const {
    if (!f1.next()) {
      // This is the simple, efficient case!
      return f1(gd, kT, f2(gd, kT, data));
    }
    // This does some extra work to save the energies of each term in
    // the sum, just in case we want to print it!
    VectorXd f2data(f2(gd, kT, data));
    VectorXd f1f2data(f1.justMe(gd, kT, f2data));
    double e = gd.dvolume*f1f2data.sum();
    f1.set_last_energy(e);
    Functional *nxt = f1.next();
    while (nxt) {
      f1f2data += nxt->justMe(gd, kT, f2data);
      double etot = gd.dvolume*f1f2data.sum();
      nxt->set_last_energy(etot - e);
      e = etot;
      nxt = nxt->next();
    }
    return f1f2data;
  }
  double transform(double kT, double n) const {
    return f1(kT, f2(kT, n));
  }
  double derive(double kT, double n) const {
    return f1.derive(kT, f2(kT, n))*f2.derive(kT, n);
  }
  double d_by_dT(double kT, double n) const {
    double f2n = f2(kT, n);
    return f1.d_by_dT(kT, f2n) + f2.d_by_dT(kT,n)*f1.derive(kT, f2n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return f2.grad(f1.grad(ingrad, f2(x), ispgrad), x, ispgrad);
  }
  Functional grad_T(const Functional &ingrad) const {
    // FIXME: The following is even less correct than below, but is
    // more efficient, as things turn out.  It assumes both that f2
    // doesn't change from real-space to some other space, and also
    // that nothing is non-local with regard to temperature.
    return f2.grad_T(ingrad)*f1.grad(Functional(1), f2, false) + f1.grad_T(ingrad)(f2);

    // FIXME: The following isn't correct if f2 changes the size of
    // the vector.  But would that even be potentially correct when
    // we're talking about temperature?
    return f1.grad(f2.grad_T(ingrad), f2, false) + f1.grad_T(ingrad)(f2);
  }
  void grad(const GridDescription &gd, double kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid outgrad1(gd);
    outgrad1.setZero();
    f1.grad(gd, kT, f2(gd, kT, data), ingrad, &outgrad1, 0);
    f2.grad(gd, kT, data, outgrad1, outgrad, outpgrad);
  }
  void print_summary(const char *prefix, double e, std::string name) const {
    f1.print_summary(prefix, e, name);
  }
  bool I_have_analytic_grad() const {
    return f1.I_have_analytic_grad() && f2.I_have_analytic_grad();
  }
private:
  Functional f1, f2;
};

Functional Functional::operator()(const Functional &f) const {
  if (f.I_am_zero() && I_give_zero_for_zero()) return f;
  return Functional(new ChainRuleType(*this, f));
}

class QuotientRuleType : public FunctionalInterface {
public:
  QuotientRuleType(const Functional &fa, const Functional &fb) : f1(fa), f2(fb) {}
  bool I_am_local() const {
    return f1.I_am_local() && f2.I_am_local();
  }
  bool I_preserve_homogeneous() const {
    return f1.I_preserve_homogeneous() && f2.I_preserve_homogeneous();
  }
  bool I_am_constant_wrt_x() const {
    return f1.I_am_constant_wrt_x() && f2.I_am_constant_wrt_x();
  }

  VectorXd transform(const GridDescription &gd, double kT, const VectorXd &data) const {
    return f1(gd, kT, data).cwise()/f2(gd, kT, data);
  }
  double transform(double kT, double n) const {
    return f1(kT, n)/f2(kT, n);
  }
  double derive(double kT, double n) const {
    double f2n = f2(kT, n);
    return f1.derive(kT, n)/f2n - f1(kT, n)*f2.derive(kT, n)/f2n/f2n;
  }
  double d_by_dT(double kT, double n) const {
    double f2n = f2(kT, n);
    return f1.d_by_dT(kT, n)/f2n - f1(kT, n)*f2.d_by_dT(kT, n)/(f2n*f2n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return f1.grad(ingrad/f2(x), x, ispgrad) - f2.grad(f1(x)*ingrad/sqr(f2(x)), x, ispgrad);
  }
  Functional grad_T(const Functional &ingrad) const {
    return f1.grad_T(ingrad)/f2 - f1/sqr(f2)*f2.grad_T(ingrad);
  }
  void grad(const GridDescription &gd, double kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    VectorXd out2 = f2(gd, kT, data);
    f1.grad(gd, kT, data, ingrad.cwise()/out2, outgrad, outpgrad);
    f2.grad(gd, kT, data, (ingrad.cwise()*f1(gd, kT, data)).cwise()/((-out2).cwise()*out2),
            outgrad, outpgrad);
  }
  bool I_have_analytic_grad() const {
    return f1.I_have_analytic_grad() && f2.I_have_analytic_grad();
  }
private:
  Functional f1, f2;
};

Functional Functional::operator/(const Functional &f) const {
  if (I_am_zero()) return *this;
  if (f.I_am_one()) return *this;
  return Functional(new QuotientRuleType(*this, f));
}

class ProductRuleType : public FunctionalInterface {
public:
  ProductRuleType(const Functional &fa, const Functional &fb) : f1(fa), f2(fb) {}
  bool I_am_local() const {
    return f1.I_am_local() && f2.I_am_local();
  }
  bool I_preserve_homogeneous() const {
    return f1.I_preserve_homogeneous() && f2.I_preserve_homogeneous();
  }
  bool I_am_constant_wrt_x() const {
    return f1.I_am_constant_wrt_x() && f2.I_am_constant_wrt_x();
  }

  VectorXd transform(const GridDescription &gd, double kT, const VectorXd &data) const {
    return f1(gd, kT, data).cwise()*f2(gd, kT, data);
  }
  double transform(double kT, double n) const {
    return f1(kT, n)*f2(kT, n);
  }
  double derive(double kT, double n) const {
    return f1(kT, n)*f2.derive(kT, n) + f1.derive(kT, n)*f2(kT, n);
  }
  double d_by_dT(double kT, double n) const {
    return f1(kT, n)*f2.d_by_dT(kT, n) + f1.d_by_dT(kT, n)*f2(kT, n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return f2.grad(f1(x)*ingrad, x, ispgrad) + f1.grad(f2(x)*ingrad, x, ispgrad);
  }
  Functional grad_T(const Functional &ingrad) const {
    return f1.grad_T(ingrad)*f2 - f1*f2.grad_T(ingrad);
  }
  void grad(const GridDescription &gd, double kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    f1.grad(gd, kT, data, ingrad.cwise()*f2(gd, kT, data), outgrad, outpgrad);
    f2.grad(gd, kT, data, ingrad.cwise()*f1(gd, kT, data), outgrad, outpgrad);
  }
  bool I_have_analytic_grad() const {
    return f1.I_have_analytic_grad() && f2.I_have_analytic_grad();
  }
private:
  Functional f1, f2;
};

Functional Functional::operator*(const Functional &f) const {
  if (I_am_one()) return f;
  if (f.I_am_one()) return *this;
  if (I_am_zero()) return *this;
  if (f.I_am_zero()) return f;
  return Functional(new ProductRuleType(*this, f));
}


class LogType : public FunctionalInterface {
public:
  LogType() {}

  VectorXd transform(const GridDescription &, double, const VectorXd &data) const {
    return data.cwise().log();
  }
  double transform(double, double n) const {
    return log(n);
  }
  double derive(double, double n) const {
    return 1/n;
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool) const {
    return ingrad/x;
  }
  Functional grad_T(const Functional &) const {
    return Functional(0.0);
  }
  void grad(const GridDescription &, double, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += ingrad.cwise()/data;
    if (outpgrad) *outpgrad += ingrad.cwise()/data;
  }
};

Functional log(const Functional &f) {
  return Functional(new LogType())(f);
}


class ExpType : public FunctionalInterface {
public:
  ExpType() {}

  VectorXd transform(const GridDescription &, double, const VectorXd &data) const {
    return data.cwise().exp();
  }
  double transform(double, double n) const {
    return exp(n);
  }
  double derive(double, double n) const {
    return exp(n);
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool) const {
    return ingrad*Functional(new ExpType())(x);
  }
  Functional grad_T(const Functional &) const {
    return Functional(0.0);
  }
  void grad(const GridDescription &, double, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += ingrad.cwise() * data.cwise().exp();
    if (outpgrad) *outpgrad += ingrad.cwise() * data.cwise().exp();
  }
};

Functional exp(const Functional &f) {
  return Functional(new ExpType())(f);
}


class AbsType : public FunctionalInterface {
public:
  AbsType() {}

  VectorXd transform(const GridDescription &, double, const VectorXd &data) const {
    return data.cwise().abs();
  }
  double transform(double, double n) const {
    return fabs(n);
  }
  double derive(double, double n) const {
    return n/fabs(n);
  }
  double d_by_dT(double, double) const {
    return 0;
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool) const {
    return ingrad*x/abs(x);
  }
  Functional grad_T(const Functional &) const {
    return Functional(0.0);
  }
  void grad(const GridDescription &, double, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += ingrad.cwise() * (data.cwise() / data.cwise().abs());
    if (outpgrad) *outpgrad += ingrad.cwise() * (data.cwise() / data.cwise().abs());
  }
};

Functional abs(const Functional &f) {
  return Functional(new AbsType())(f);
}


Functional Functional::operator-(const Functional &f) const {
  return *this + -f;
}

Functional Functional::operator-() const {
  return -1* *this;
}

Functional sqr(const Functional &f) {
  return Pow(2)(f);
}


class Constraint : public FunctionalInterface {
public:
  Constraint(const Grid &g, const Functional &y) : constraint(g), f(y) {};
  bool I_am_local() const {
    return f.I_am_local();
  }
  bool I_am_constant_wrt_x() const {
    return f.I_am_constant_wrt_x();
  }
  bool I_preserve_homogeneous() const {
    return f.I_preserve_homogeneous();
  }
  bool I_am_homogeneous() const {
    return f.I_am_homogeneous();
  }
  bool I_am_zero() const {
    return f.I_am_zero();
  }
  bool I_am_one() const {
    return f.I_am_one();
  }

  double integral(const GridDescription &gd, double kT, const VectorXd &x) const {
    return f.integral(gd, kT, x);
  }
  VectorXd transform(const GridDescription &gd, double kT, const VectorXd &data) const {
    return f(gd, kT, data);
  }
  double transform(double kT, double n) const {
    return f(kT, n);
  }
  double derive(double kT, double n) const {
    return f.derive(kT, n);
  }
  double d_by_dT(double kT, double n) const {
    return f.d_by_dT(kT, n);
  }
  Functional grad(const Functional &ingrad, const Functional &x, bool ispgrad) const {
    return Functional(constraint)*f.grad(ingrad, x, ispgrad);
  }
  Functional grad_T(const Functional &ingrad) const {
    return Functional(constraint)*f.grad_T(ingrad);
  }
  void pgrad(const GridDescription &gd, double kT, const VectorXd &data,
             const VectorXd &ingrad, VectorXd *outpgrad) const {
    VectorXd mypgrad(data);
    mypgrad.setZero();
    f.pgrad(gd, kT, data, ingrad, &mypgrad);
    *outpgrad += constraint.cwise() * mypgrad;
  }
  void grad(const GridDescription &gd, double kT, const VectorXd &data,
            const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
    if (outpgrad) {
      VectorXd mygrad(data), mypgrad(data);
      mygrad.setZero();
      mypgrad.setZero();
      f.grad(gd, kT, data, ingrad, &mygrad, &mypgrad);
      *outgrad += constraint.cwise() * mygrad;
      *outpgrad += constraint.cwise() * mypgrad;
    } else {
      VectorXd mygrad(data);
      mygrad.setZero();
      f.grad(gd, kT, data, ingrad, &mygrad, 0);
      *outgrad += constraint.cwise() * mygrad;
    }
  }
  void print_summary(const char *prefix, double e, std::string name) const {
    f.print_summary(prefix, e, name);
  }
  bool I_have_analytic_grad() const {
    return f.I_have_analytic_grad();
  }
private:
  const Grid constraint;
  Functional f;
};

Functional constrain(const Grid &g, Functional f) {
  return Functional(new Constraint(g, f));
}
