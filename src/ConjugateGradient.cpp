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

#include "LineMinimizer.h"
#include <stdio.h>

class ConjugateGradientType : public MinimizerInterface {
protected:
  double step, orig_step;
  VectorXd direction, oldgrad;
  LineMinimizer linmin;
  double oldgradsqr;
public:
  ConjugateGradientType(Functional f, const GridDescription &gdin, VectorXd *data, LineMinimizer lm,
                        double stepsize = 0.1)
    : MinimizerInterface(f, gdin, data), step(stepsize), orig_step(step), direction(*data), oldgrad(*data), linmin(lm) {
    direction.setZero();
    oldgrad.setZero();
    oldgradsqr = 0;
  }
  void minimize(Functional newf, const GridDescription &gdnew, VectorXd *newx = 0) {
    step = orig_step;
    MinimizerInterface::minimize(newf, gdnew, newx);
    if (newx) {
      direction = *newx;
      oldgrad = *newx;
    }
    direction.setZero();
    oldgrad.setZero();
    oldgradsqr = 0;
  }

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};

class PreconditionedConjugateGradientType : public ConjugateGradientType {
public:
  PreconditionedConjugateGradientType(Functional f, const GridDescription &gdin,
                                    VectorXd *data, LineMinimizer lm, double stepsize = 0.1)
    : ConjugateGradientType(f, gdin, data, lm, stepsize) {}

  bool improve_energy(bool verbose = false);
};

bool ConjugateGradientType::improve_energy(bool verbose) {
  iter++;
  //printf("I am running ConjugateGradient::improve_energy\n");
  const double E0 = energy();
  double gdotd;
  {
    const VectorXd g = -grad();
    // Let's immediately free the cached gradient stored internally!
    invalidate_cache();

    // Note: my notation vaguely follows that of
    // [wikipedia](http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method).
    // I use the Polak-Ribiere method, with automatic direction reset.
    // Note that we could save some memory by using Fletcher-Reeves, and
    // it seems worth implementing that as an option for
    // memory-constrained problems (then we wouldn't need to store oldgrad).
    double beta = g.dot(g - oldgrad)/oldgradsqr;
    oldgrad = g;
    if (beta < 0 || beta != beta || oldgradsqr == 0) beta = 0;
    oldgradsqr = oldgrad.dot(oldgrad);
    direction = g + beta*direction;
    gdotd = oldgrad.dot(direction);
    if (gdotd < 0) {
      direction = oldgrad; // If our direction is uphill, reset to gradient.
      if (verbose) printf("reset to gradient...\n");
      gdotd = oldgrad.dot(direction);
    }
  }

  Minimizer lm = linmin(f, gd, x, direction, -gdotd, &step);
  for (int i=0; i<100 && lm.improve_energy(verbose); i++) {
    if (verbose) lm.print_info("\t");
  }
  if (verbose) {
    //lm->print_info();
    print_info();
    printf("grad*oldgrad = %g\n", grad().dot(direction)/gdotd);
  }
  return (energy() < E0);
}

void ConjugateGradientType::print_info(const char *prefix) const {
  MinimizerInterface::print_info(prefix);
  printf("%sstep = %g\n", prefix, step);
}

bool PreconditionedConjugateGradientType::improve_energy(bool verbose) {
  iter++;
  //printf("I am running ConjugateGradient::improve_energy\n");
  const double E0 = energy();
  double beta;
  {
    // Note: my notation vaguely follows that of
    // [wikipedia](http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method).
    // I use the Polak-Ribiere method, with automatic direction reset.
    // Note that we could save some memory by using Fletcher-Reeves, and
    // it seems worth implementing that as an option for
    // memory-constrained problems (then we wouldn't need to store oldgrad).
    pgrad(); // compute pgrad first, since that computes both.
    beta = -pgrad().dot(-grad() - oldgrad)/oldgradsqr;
    oldgrad = -grad();
    if (beta < 0 || beta != beta || oldgradsqr == 0) beta = 0;
    if (verbose) printf("beta = %g\n", beta);
    oldgradsqr = -pgrad().dot(oldgrad);
    direction = -pgrad() + beta*direction;
    // Let's immediately free the cached gradient stored internally!
    invalidate_cache();
  } // free g and pg!

  const double gdotd = oldgrad.dot(direction);

  Minimizer lm = linmin(f, gd, x, direction, -gdotd, &step);
  for (int i=0; i<100 && lm.improve_energy(verbose); i++) {
    if (verbose) lm.print_info("\t");
  }
  if (verbose) {
    //lm->print_info();
    print_info();
    printf("grad*oldgrad = %g\n", grad().dot(direction)/gdotd);
  }
  return (energy() < E0 || beta != 0);
}


Minimizer ConjugateGradient(Functional f, const GridDescription &gdin, VectorXd *data,
                          LineMinimizer lm, double stepsize) {
  return Minimizer(new ConjugateGradientType(f, gdin, data, lm, stepsize));
}

Minimizer PreconditionedConjugateGradient(Functional f, const GridDescription &gdin, VectorXd *data,
                                        LineMinimizer lm, double stepsize) {
  return Minimizer(new PreconditionedConjugateGradientType(f, gdin, data, lm, stepsize));
}
