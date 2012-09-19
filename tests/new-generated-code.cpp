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

#include "new/Functional.h"
#include "new-generated-haskell/volume_minus_one_sqr.h"
#include "new-generated-haskell/integrate_sqr.h"
#include "new-generated-haskell/WhiteBear.h"
#include "handymath.h"

int errors = 0;

double a = 5;

class volminusonesqr : public Functional {
public:
  volminusonesqr() {}
  double energy(const Vector &x) const {
    Vector a1 = x.slice(0,3);
    Vector a2 = x.slice(3,3);
    Vector a3 = x.slice(6,3);
    V = 0;
    for (int i=0;i<3;i++) {
      V += a1[i]*(a2[(i+1)%3]*a3[(i+2)%3] - a2[(i+2)%3]*a3[(i+1)%3]);
    }
    return (V-1)*(V-1);
  }
  double energy_per_volume(const Vector &x) const {
    return 0; // What does this mean?
  }
  double denergy_per_volume_dx(const Vector &x) const {
    return 0; // What does this mean?
  }
  Vector grad(const Vector &x) const {
    Vector out(9);
    for (int i=0;i<9;i++) out[i] = 0;
    return out;
  }
  void printme(const char *prefix) const {
    printf("%s V = %g\n", prefix, V);
  }
private:
	mutable double V;
};

class nsqr : public Functional {
public:
  nsqr() {}
  double energy(const Vector &x) const {
    double Nx = x[0];
    double Ny = x[1];
    double Nz = x[2];
    double a1 = x[3];
    double a2 = x[4];
    double a3 = x[5];
    printf("foobar foobar\n");
    Vector n = x.slice(6,Nx*Ny*Nz);
    double volume = a1*a2*a3; // only works for simple cubic lattice
    double dV = volume/Nx/Ny/Nz;
    double E = 0;
    for (int i=0;i<Nx*Ny*Nz;i++) {
      E += dV*n[i]*n[i];
    }
    return E;
  }
  double energy_per_volume(const Vector &x) const {
    return 0; // What does this mean?
  }
  double denergy_per_volume_dx(const Vector &x) const {
    return 0; // What does this mean?
  }
  Vector grad(const Vector &x) const {
    double Nx = x[0];
    double Ny = x[1];
    double Nz = x[2];
    double a1 = x[3];
    double a2 = x[4];
    double a3 = x[5];
    printf("foobar bazbar\n");
    Vector n = x.slice(6,Nx*Ny*Nz);
    double volume = a1*a2*a3; // only works for simple cubic lattice
    double dV = volume/Nx/Ny/Nz;
    double E = 0;
    Vector g(x.get_size());
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
                         const Functional &f1, const Functional &f2, Vector v,
                         double fraccuracy = 1e-15) {
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  printf("First energy:\n");
  double f1n = f1.energy(v);
  print_double("first energy is:               ", f1n);
  printf("\n");
  f1.printme("  ");
  printf("Second energy:\n");
  double f2n = f2.energy(v);
  print_double("second energy is:              ", f2n);
  printf("\n");
  f2.printme("  ");
  if (fabs(f1n/f2n - 1) > fraccuracy) {
    printf("E1 = %g\n", f1n);
    printf("E2 = %g\n", f2n);
    printf("FAIL: Error in f(n) is %g\n", f1n/f2n - 1);
    errors++;
  }
  Vector g1 = f1.grad(v);
  Vector g2 = f2.grad(v);
  double err = (g1-g2).norm();
  double mag = g1.norm();
  if (err/mag > fraccuracy) {
    printf("FAIL: Error in grad %s is %g as a fraction of %g\n", name, err/mag, mag);
    errors++;
  }
}

int main(int, char **argv) {

  volminusonesqr vm1s;
  volume_minus_one_sqr vm1sgen;
  Vector lat(9);
  lat[0] = 5; lat[1] = 0; lat[2] = 0;
  lat[3] = 0; lat[4] = 5; lat[5] = 0;
  lat[6] = 0; lat[7] = 0; lat[8] = 5;
  compare_functionals("volume minus 1 squared", vm1s, vm1sgen, lat, 1e-10);

  Vector n(1000);
  for (int i=0;i<1000;i++) n[i] = 0.5;
  Vector x = integrate_sqr().createInput(10,10,10,
                                         3, 3, 3,
                                         n);
  compare_functionals("integrate_sqr", nsqr(), integrate_sqr(), x, 1e-10);

  if (errors == 0) printf("\n%s passes!\n", argv[0]);
  else printf("\n%s fails %d tests!\n", argv[0], errors);
  return errors;
}
