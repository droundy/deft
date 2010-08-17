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

#include "Minimizer.h"
#include <stdio.h>

class DownhillType : public MinimizerInterface {
protected:
  double nu, orig_nu;
public:
  DownhillType(Functional f, const GridDescription &gdin, VectorXd *data, double viscosity=0.1)
    : MinimizerInterface(f, gdin, data), nu(viscosity), orig_nu(viscosity) {}
  ~DownhillType() {}
  void minimize(Functional newf, const GridDescription &gdnew, VectorXd *newx = 0) {
    nu = orig_nu;
    MinimizerInterface::minimize(newf, gdnew, newx);
  }

  bool improve_energy(bool verbose = false);
  void print_info(const char *prefix="") const;
};

class PreconditionedDownhillType : public DownhillType {
public:
  PreconditionedDownhillType(Functional f, const GridDescription &gdin, VectorXd *data, double viscosity=0.1)
    : DownhillType(f, gdin, data, viscosity) {}
  ~PreconditionedDownhillType() {}

  bool improve_energy(bool verbose = false);
};

bool DownhillType::improve_energy(bool verbose) {
  iter++;
  const VectorXd g = grad();
  // Let's immediately free the cached gradient stored internally!
  invalidate_cache();
  // We waste some memory storing oldx, but avoids roundoff weirdness
  // of trying to add nu*g back to *x, which won't always get us back
  // to the same value.
  const VectorXd oldx = *x;
  double old_energy = energy();
  *x -= nu*g;
  invalidate_cache(); // Must always remember this!!!
  int num_tries = 0;
  while (energy() > old_energy) {
    nu *= 0.5;
    *x = oldx - nu*g;
    invalidate_cache();
    if (num_tries++ > 30) {
      printf("Downhill giving up after %d tries...\n", num_tries);
      return false; // It looks like we can't do any better with this algorithm.
    }
  }
  nu *= 1.1;
  if (verbose) {
    //lm->print_info();
    print_info();
  }
  return true;
}

void DownhillType::print_info(const char *prefix) const {
  MinimizerInterface::print_info(prefix);
  printf("%snu = %g\n", prefix, nu);
}

bool PreconditionedDownhillType::improve_energy(bool verbose) {
  iter++;
  const VectorXd g = pgrad();
  // Let's immediately free the cached gradient stored internally!
  invalidate_cache();
  // We waste some memory storing oldx, but avoids roundoff weirdness
  // of trying to add nu*g back to *x, which won't always get us back
  // to the same value.
  const VectorXd oldx = *x;
  double old_energy = energy();
  *x -= nu*g;
  invalidate_cache(); // Must always remember this!!!
  int num_tries = 0;
  while (energy() > old_energy) {
    nu *= 0.5;
    *x = oldx - nu*g;
    invalidate_cache();
    if (num_tries++ > 30) {
      printf("PreconditionedDownhill giving up after %d tries, with gradient %g...\n", num_tries, g.norm());
      return false; // It looks like we can't do any better with this algorithm.
    }
  }
  nu *= 1.1;
  if (verbose) {
    //lm->print_info(iter);
    print_info();
  }
  return true;
}

Minimizer Downhill(Functional f, const GridDescription &gdin, VectorXd *data, double viscosity) {
  return Minimizer(new DownhillType(f, gdin, data, viscosity));
}

Minimizer PreconditionedDownhill(Functional f, const GridDescription &gdin, VectorXd *data, double viscosity) {
  return Minimizer(new PreconditionedDownhillType(f, gdin, data, viscosity));
}
