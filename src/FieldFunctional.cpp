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

void FieldFunctionalInterface::print_summary(const char *prefix, double e, const char *name) const {
  if (name) printf("%s%25s =", prefix, name);
  else printf("%s%25s =", prefix, "UNKNOWN");
  print_double("", e);
  printf("\n");
}

bool FieldFunctional::run_finite_difference_test(const char *testname, const Grid &x,
                                                 const VectorXd *direction) const {
  printf("\nRunning finite difference test on %s:\n", testname);
  const GridDescription gd(x.description());

  VectorXd my_grad(x);
  double Eold = integral(x);
  my_grad.setZero();
  integralgrad(x, &my_grad);
  VectorXd my_direction(my_grad);
  if (direction) my_direction = *direction;
  my_direction /= my_direction.norm();

  const double lderiv = my_direction.dot(my_grad);
  if (lderiv == 0.0) {
    printf("Gradient is zero...\n");
    return true;
  }
  {
    // compare the gradient computed by the two functions grad() and
    // grad_and_pgrad()
    VectorXd my_pgrad(x);
    my_grad.setZero();
    my_pgrad.setZero();
    integralgrad(x, &my_grad, &my_pgrad);
    const double lderiv_new = my_direction.dot(my_grad);
    if (fabs(lderiv_new/lderiv - 1) > 1e-12) {
      printf("\n*** WARNING!!! INCONSISTENT GRADIENTS! ***\n");
      printf("Different gradient in energy_and_grad() and grad_and_pgrad()\n");
      printf("Fractional error is %g\n\n", lderiv_new/lderiv - 1);
      return false;
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
    const double Eplus= integral(gd, x + epsilon*my_direction);
    const double Eminus=integral(gd, x - epsilon*my_direction);

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
    return false;
  }
  if (best_ratio_error < 1e-10) {
    printf("Passed on basis of a gradient accuracy of (%g) (with scaling %g).\n",
           best_ratio_error, min);
    return false;
  }
  if (min < 1e-5 && best_ratio_error < 0.01) {
    printf("Passed on basis of seriously nice scaling %g (and low accuracy %g).\n",
           min, best_ratio_error);
    return false;
  }
  printf("FAIL: Failed with scaling ratio of %g and gradient accuracy of %g\n",
         min, best_ratio_error);
  return true;
}

double FieldFunctional::print_iteration(const char *prefix, int iter) const {
  printf("%s==============\n", prefix);
  printf("%sIteration %4d\n", prefix, iter);
  printf("%s==============\n", prefix);

  print_summary(prefix, itsCounter->last_energy, get_name());
  const FieldFunctional *nxt = this;
  double etot = 0;
  while (nxt) {
    etot += nxt->itsCounter->last_energy;
    //printf("In print_iteration summarizing %p\n", nxt);
    //nxt->print_summary(prefix, nxt->itsCounter->last_energy);
    nxt = nxt->next();
  }
  printf("%s%25s =", prefix, "total energy");
  print_double("", etot);
  printf("\n");
  return etot;
}

void FieldFunctional::print_summary(const char *prefix, double energy, const char *name) const {
  if (get_name()) name = get_name();
  if (!next()) {
    itsCounter->ptr->print_summary(prefix, energy, name);
    return;
  }
  const FieldFunctional *nxt = this;
  double etot = 0;
  while (nxt) {
    if (nxt->get_name()) name = nxt->get_name();
    etot += nxt->itsCounter->last_energy;
    nxt->itsCounter->ptr->print_summary(prefix, nxt->itsCounter->last_energy, name);
    nxt = nxt->next();
  }
  //assert(fabs(etot - energy) < 1e-6);
}

class IdentityType : public FieldFunctionalInterface {
public:
  IdentityType() {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    return data;
  }
  double transform(double n) const {
    return n;
  }
  double grad(double) const {
    return 1;
  }

  void grad(const GridDescription &, const VectorXd &, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    *outgrad += ingrad;
    if (outpgrad) *outpgrad += ingrad; // FIXME: propogate preconditioning
  }
};

FieldFunctional Identity() { return FieldFunctional(new IdentityType()); }

class Constant : public FieldFunctionalInterface {
public:
  Constant(double x) : c(x) {}

  VectorXd transform(const GridDescription &, const VectorXd &data) const {
    return c*VectorXd::Ones(data.rows());
  }
  double transform(double) const {
    return c;
  }
  double grad(double) const {
    return 0;
  }

  void grad(const GridDescription &, const VectorXd &, const VectorXd &, VectorXd *, VectorXd *) const {
  }
private:
  double c;
};

FieldFunctional::FieldFunctional(double x) : itsCounter(0) {
  init(new Constant(x), 0);
}

class ConstantField : public FieldFunctionalInterface {
public:
  ConstantField(const VectorXd &x) : c(x) {}

  VectorXd transform(const GridDescription &, const VectorXd &) const {
    return c;
  }
  double transform(double) const {
    return c.sum()/c.rows();
  }
  double grad(double) const {
    return 0;
  }

  void grad(const GridDescription &, const VectorXd &, const VectorXd &, VectorXd *, VectorXd *) const {
  }
private:
  VectorXd c;
};

FieldFunctional::FieldFunctional(const VectorXd &x) : itsCounter(0) {
  init(new ConstantField(x), 0);
}

class ChainRuleType : public FieldFunctionalInterface {
public:
  ChainRuleType(const FieldFunctional &fa, const FieldFunctional &fb) : f1(fa), f2(fb) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    if (!f1.next()) {
      // This is the simple, efficient case!
      return f1(gd, f2(gd, data));
    }
    // This does some extra work to save the energies of each term in
    // the sum, just in case we want to print it!
    VectorXd f2data(f2(gd,data));
    VectorXd f1f2data(f1.justMe(gd, f2data));
    double e = gd.dvolume*f1f2data.sum();
    f1.set_last_energy(e);
    FieldFunctional *nxt = f1.next();
    while (nxt) {
      f1f2data += nxt->justMe(gd, f2data);
      double etot = gd.dvolume*f1f2data.sum();
      nxt->set_last_energy(etot - e);
      e = etot;
      nxt = nxt->next();
    }
    return f1f2data;
  }
  double transform(double n) const {
    return f1(f2(n));
  }
  double grad(double n) const {
    return f1.grad(f2(n))*f2.grad(n);
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    Grid outgrad1(gd);
    outgrad1.setZero();
    f1.grad(gd, f2(gd, data), ingrad, &outgrad1, 0);
    f2.grad(gd, data, outgrad1, outgrad, outpgrad);
  }
  void print_summary(const char *prefix, double e, const char *name) const {
    f1.print_summary(prefix, e, name);
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
  double grad(double n) const {
    return f1.grad(n)/f2(n) - f1(n)*f2.grad(n)/f2(n)/f2(n);
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
  double grad(double n) const {
    return f1(n)*f2.grad(n) + f1.grad(n)*f2(n);
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


class ScalarMinusType : public FieldFunctionalInterface {
public:
  ScalarMinusType(const FieldFunctional &fa, double xx) : f(fa), x(xx) {}

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return x*VectorXd::Ones(gd.NxNyNz) - f(gd, data);
  }
  double transform(double n) const {
    return x - f(n);
  }
  double grad(double n) const {
    return -f.grad(n);
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
  double grad(double n) const {
    return f.grad(n)/f(n);
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
  double grad(double n) const {
    return 2*f(n)*f.grad(n);
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


class Constraint : public FieldFunctionalInterface {
public:
  Constraint(const Grid &g, const FieldFunctional &y) : constraint(g), f(y) {};

  VectorXd transform(const GridDescription &gd, const VectorXd &data) const {
    return f(gd, data);
  }
  double transform(double n) const {
    return f(n);
  }
  double grad(double n) const {
    return f.grad(n);
  }

  void grad(const GridDescription &gd, const VectorXd &data, const VectorXd &ingrad,
            VectorXd *outgrad, VectorXd *outpgrad) const {
    if (outpgrad) {
      VectorXd mygrad(data), mypgrad(data);
      mygrad.setZero();
      mypgrad.setZero();
      f.grad(gd, data, ingrad, &mygrad, &mypgrad);
      *outgrad += constraint.cwise() * mygrad;
      *outpgrad += constraint.cwise() * mypgrad;
    } else {
      VectorXd mygrad(data);
      mygrad.setZero();
      f.grad(gd, data, ingrad, &mygrad, 0);
      *outgrad += constraint.cwise() * mygrad;
    }
  }
private:
  const Grid constraint;
  const FieldFunctional f;
};

FieldFunctional constrain(const Grid &g, FieldFunctional f) {
  return FieldFunctional(new Constraint(g, f));
}
