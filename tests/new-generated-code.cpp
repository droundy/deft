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

#include "new/NewFunctional.h"
#include "new-generated-haskell/volume_minus_one_sqr.h"
#include "new-generated-haskell/integrate_sqr.h"
#include "new-generated-haskell/WhiteBear.h"
#include "handymath.h"

int errors = 0;

class nsqr : public NewFunctional {
public:
  nsqr() {}
  double true_energy() const {
    double Nx = data[0];
    double Ny = data[1];
    double Nz = data[2];
    double a1 = data[3];
    double a2 = data[4];
    double a3 = data[5];
    Vector n = data.slice(6,Nx*Ny*Nz);
    double volume = a1*a2*a3; // only works for simple orthorhombic lattice
    double dV = volume/Nx/Ny/Nz;
    double E = 0;
    for (int i=0;i<Nx*Ny*Nz;i++) {
      E += dV*n[i]*n[i];
    }
    return E;
  }
  Vector grad() const {
    double Nx = data[0];
    double Ny = data[1];
    double Nz = data[2];
    double a1 = data[3];
    double a2 = data[4];
    double a3 = data[5];
    Vector n = data.slice(6,Nx*Ny*Nz);
    double volume = a1*a2*a3; // only works for simple cubic lattice
    double dV = volume/Nx/Ny/Nz;
    double E = 0;
    Vector g(data.get_size());
    g *= 0.0;
    for (int i=0;i<Nx*Ny*Nz;i++) {
      g[i+6] += 2*dV*n[i];
    }
    return g;
  }
  void printme(const char *prefix) const {
  }
};

void compare_functionals(const char *name,
                         const NewFunctional &f1, const NewFunctional &f2,
                         double fraccuracy = 1e-15) {
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  printf("First energy:\n");
  double f1n = f1.energy();
  print_double("first energy is:               ", f1n);
  printf("\n");
  f1.printme("  ");
  printf("Second energy:\n");
  double f2n = f2.energy();
  print_double("second energy is:              ", f2n);
  printf("\n");
  f2.printme("  ");
  if (fabs(f1n/f2n - 1) > fraccuracy) {
    printf("E1 = %g\n", f1n);
    printf("E2 = %g\n", f2n);
    printf("FAIL: Error in f(n) is %g\n", f1n/f2n - 1);
    errors++;
  }
  Vector g1 = f1.grad();
  Vector g2 = f2.grad();
  double err = (g1-g2).norm();
  double mag = g1.norm();
  if (err/mag > fraccuracy) {
    printf("FAIL: Error in grad %s is %g as a fraction of %g\n", name, err/mag, mag);
    errors++;
  }
}

int main(int, char **argv) {

  Vector n(1000);
  n = 0.5;
  integrate_sqr haskell_generated;
  haskell_generated.alloc(10,10,10);
  haskell_generated.a1() = 5;
  haskell_generated.a2() = 5;
  haskell_generated.a3() = 5;
  haskell_generated.nn() = n;

  nsqr hand_coded;
  hand_coded.copy_input_data_for_test(haskell_generated); // use same input for both
  compare_functionals("integrate_sqr", hand_coded, haskell_generated, 1e-10);

  if (errors == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], errors);
  return errors;
}
