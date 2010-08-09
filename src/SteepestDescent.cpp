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

#include "SteepestDescent.h"
#include <stdio.h>

bool SteepestDescent::improve_energy(bool verbose) {
  iter++;
  //printf("I am running SteepestDescent::improve_energy\n");
  const double E0 = energy();
  const VectorXd d = -grad();
  // Let's immediately free the cached gradient stored internally!
  invalidate_cache();

  Minimizer *lm = linmin(f, x, d, -d.dot(d), &step);
  for (int i=0; i<100 && lm->improve_energy(verbose); i++) {
    if (verbose) lm->print_info("\t");
  }
  if (verbose) {
    //lm->print_info();
    print_info();
  }
  delete lm;
  return (energy() < E0);
}

void SteepestDescent::print_info(const char *prefix) const {
  Minimizer::print_info(prefix);
  printf("%sstep = %g\n", prefix, step);
}

bool PreconditionedSteepestDescent::improve_energy(bool verbose) {
  iter++;
  //printf("I am running PreconditionedSteepestDescent::improve_energy\n");
  const double E0 = energy();
  const VectorXd d = -pgrad();
  // Let's immediately free the cached gradient stored internally!
  invalidate_cache();

  Minimizer *lm = linmin(f, x, d, d.dot(grad()), &step);
  for (int i=0; i<100 && lm->improve_energy(verbose); i++) {
    if (verbose) lm->print_info("\t");
  }
  if (verbose) {
    //lm->print_info();
    print_info();
  }
  return (energy() < E0);
}

void PreconditionedSteepestDescent::print_info(const char *prefix) const {
  Minimizer::print_info(prefix);
  printf("%sstep = %g\n", prefix, step);
}
