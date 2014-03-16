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

#include <stdio.h>
#include <time.h>
#include "new/QuadraticFast.h"
#include "new/QuadraticN0Fast.h"
#include "new/QuadraticGaussianFast.h"
#include "new/LogN0Fast.h"
#include "new/Phi1Fast.h"
#include "new/Phi2Fast.h"
#include "new/Phi3Fast.h"
#include "new/SPhi1Fast.h"
#include "new/SPhi2Fast.h"
#include "new/SPhi3Fast.h"
#include "new/SFMTFluidFast.h"
#include "new/HomogeneousSFMTFluidFast.h"
#include "new/ExternalPotentialTestFast.h"
#include "handymath.h"

int check_functional_value(const char *name,
                            const NewFunctional &f,
                            double energy,
                            double fraccuracy = 1e-15) {
  int errors = 0;
  printf("\n************");
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("\n* Testing %s *\n", name);
  for (unsigned i=0;i<strlen(name);i++) printf("*");
  printf("************\n\n");

  printf("findiing energy\n");
  double fv = f.energy();
  printf("found energy\n");
  print_double("Energy of Vector:  ", fv);
  printf("\n");
  f.printme("  ");

  if (!(fabs(fv/energy - 1) < fraccuracy)) {
    printf("fv = %.15g\n", fv);
    printf("expected = %.15g\n", energy);
    printf("FAIL: Error in f(n) is %g\n", fv/energy - 1);
    errors++;
  }
  return errors;
}

int main(int, char **argv) {
  int retval = 0;
  const int run_slow_tests = 1;
  //const int run_failing_tests = 1;

  int Nx = 100;
  double a = 2;

  if (run_slow_tests) {
    Quadratic q(Nx, Nx, Nx);
    q.a1() = a;
    q.a2() = a;
    q.a3() = a;
    q.konst() = 1;
    q.meanval() = 1;
    q.x() = 1.5;
    double energy = 0.5*q.konst()*a*a*a*sqr(q.x()[0]-q.meanval());
    retval += check_functional_value("Quadratic", q, energy, 2e-11);
    retval += q.run_finite_difference_test("Quadratic");
  }
  if (run_slow_tests) {
    QuadraticN0 f(Nx, Nx, Nx);
    f.a1() = a;
    f.R() = 1;
    f.a2() = a;
    f.a3() = a;
    f.konst() = 1;
    f.meanval() = 1;
    f.x() = 1.5;
    double energy = 0.5*f.konst()*a*a*a*sqr(f.x()[0]-f.meanval());
    retval += check_functional_value("Quadratic n0", f, energy, 2e-11);
    retval += f.run_finite_difference_test("Quadratic n0");
    Vector z = f.get_rz();
    for (int i=0; i<z.get_size(); i++) {
      if (z[i] > 0) {
        f.x()[i] *= 0.0001;
      }
    }
    retval += f.run_finite_difference_test("Quadratic n0 inhomogeneous");
  }
  if (run_slow_tests) {
    QuadraticGaussian f(Nx, Nx, Nx);
    f.a1() = a;
    f.a2() = a;
    f.a3() = a;
    f.sigma() = 1.0;
    f.konst() = 1;
    f.meanval() = 1;
    f.x() = 1.5;
    double energy = 0.5*f.konst()*a*a*a*sqr(f.x()[0]-f.meanval());
    retval += check_functional_value("Quadratic gaussian", f, energy, 2e-11);
    retval += f.run_finite_difference_test("Quadratic gaussian");
    Vector z = f.get_rz();
    for (int i=0; i<z.get_size(); i++) {
      if (z[i] > 0) {
        f.x()[i] *= 0.0001;
      }
    }
    retval += f.run_finite_difference_test("Quadratic gaussian inhomogeneous");
  }
  if (run_slow_tests) {
    LogN0 L0(Nx, Nx, Nx);
    L0.R() = 1;
    L0.a1() = a;
    L0.a2() = a;
    L0.a3() = a;
    L0.x() = 1.5;
    double energy = a*a*a*log(L0.x()[0]);
    retval += check_functional_value("Log n0", L0, energy, 2e-11);
    retval += L0.run_finite_difference_test("Log n0");
    for (int i=0; i<Nx*Nx*Nx/2; i++) L0.x()[i] *= 0.01;
    retval += L0.run_finite_difference_test("Log n0 inhomogeneous");
  }
  if (run_slow_tests || 1) {
    SPhi1 f(Nx, Nx, Nx);
    f.R() = 1;
    f.kT() = 1;
    f.V0() = 1;
    f.a1() = a;
    f.a2() = a;
    f.a3() = a;
    f.x() = 0.1;
    //double energy = a*a*a*log(f.x()[0]);
    //retval += check_functional_value("SPhi1", f, energy, 2e-11);
    retval += f.run_finite_difference_test("SPhi1");
    Vector z = f.get_rz();
    for (int i=0; i<z.get_size(); i++) {
      if (z[i] > 0) {
        f.x()[i] *= 0.0001;
      }
    }
    retval += f.run_finite_difference_test("SPhi1 inhomogeneous");
  }
  if (run_slow_tests) {
    SPhi2 sphi2(Nx, Nx, Nx);
    sphi2.R() = 1;
    sphi2.kT() = 1;
    sphi2.V0() = 1;
    sphi2.a1() = a;
    sphi2.a2() = a;
    sphi2.a3() = a;
    sphi2.x() = 0.1;
    //double energy = a*a*a*log(sphi2.x()[0]);
    //retval += check_functional_value("SPhi2", sphi2, energy, 2e-11);
    retval += sphi2.run_finite_difference_test("SPhi2");
    for (int i=0; i<Nx*Nx*Nx/2; i++) sphi2.x()[i] *= 0.1;
    retval += sphi2.run_finite_difference_test("SPhi2 inhomogeneous");
  }
  if (run_slow_tests) {
    SPhi3 sphi3(Nx, Nx, Nx);
    sphi3.R() = 1;
    sphi3.kT() = 1;
    sphi3.V0() = 1;
    sphi3.a1() = a;
    sphi3.a2() = a;
    sphi3.a3() = a;
    sphi3.x() = 0.1;
    //double energy = a*a*a*log(sphi3.x()[0]);
    //retval += check_functional_value("SPhi3", sphi3, energy, 2e-11);
    retval += sphi3.run_finite_difference_test("SPhi3");
    for (int i=0; i<Nx*Nx*Nx/2; i++) sphi3.x()[i] *= 0.1;
    retval += sphi3.run_finite_difference_test("SPhi3 inhomogeneous");
  }
  // And testing bits of fundamental measure theory:
  if (run_slow_tests) {
    Phi1 phi1(Nx, Nx, Nx);
    phi1.R() = 1;
    phi1.kT() = 1;
    phi1.a1() = a;
    phi1.a2() = a;
    phi1.a3() = a;
    phi1.x() = 0.1;
    //double energy = a*a*a*log(phi1.x()[0]);
    //retval += check_functional_value("Phi1", phi1, energy, 2e-11);
    retval += phi1.run_finite_difference_test("Phi1");
    for (int i=0; i<Nx*Nx*Nx/2; i++) phi1.x()[i] *= 0.1;
    retval += phi1.run_finite_difference_test("Phi1 inhomogeneous");
  }
  if (run_slow_tests) {
    Phi2 phi2(Nx, Nx, Nx);
    phi2.R() = 1;
    phi2.kT() = 1;
    phi2.a1() = a;
    phi2.a2() = a;
    phi2.a3() = a;
    phi2.x() = 0.1;
    //double energy = a*a*a*log(phi2.x()[0]);
    //retval += check_functional_value("Phi2", phi2, energy, 2e-11);
    retval += phi2.run_finite_difference_test("Phi2");
    for (int i=0; i<Nx*Nx*Nx/2; i++) phi2.x()[i] *= 0.1;
    retval += phi2.run_finite_difference_test("Phi2 inhomogeneous");
  }
  if (run_slow_tests) {
    Phi3 phi3(Nx, Nx, Nx);
    phi3.R() = 1;
    phi3.kT() = 1;
    phi3.a1() = a;
    phi3.a2() = a;
    phi3.a3() = a;
    phi3.x() = 0.1;
    //double energy = a*a*a*log(phi3.x()[0]);
    //retval += check_functional_value("Phi3", phi3, energy, 2e-11);
    retval += phi3.run_finite_difference_test("Phi3");
    for (int i=0; i<Nx*Nx*Nx/2; i++) phi3.x()[i] *= 0.1;
    retval += phi3.run_finite_difference_test("Phi3 inhomogeneous");
  }
  if (run_slow_tests) {
    ExternalPotentialTest ep(Nx, Nx, Nx);
    ep.n() = 1;
    ep.V() = 1;
    ep.a1() = a;
    ep.a2() = a;
    ep.a3() = a;
    double energy = a*a*a;
    retval += check_functional_value("ExternalPotentialTest", ep, energy, 2e-11);
    retval += ep.run_finite_difference_test("ExternalPotentialTest");
    for (int i=0; i<Nx*Nx*Nx/2; i++) {
      ep.n()[i] *= 0.1;
      ep.V()[i] = i;
    }
    retval += ep.run_finite_difference_test("ExternalPotentialTest inhomogeneous");
  }
  if (run_slow_tests) {
    HomogeneousSFMTFluid hf;
    hf.R() = 1;
    hf.V0() = 1;
    hf.kT() = 0.001;
    hf.n() = 0.3/(4*M_PI/3);
    hf.mu() = 0;
    hf.mu() = hf.d_by_dn(); // set mu based on derivative of hf

    SFMTFluid f(2, 2, Nx);
    f.a1() = a*2/Nx;
    f.a2() = a*2/Nx;
    f.a3() = a;
    f.R() = hf.R();
    f.V0() = hf.V0();
    f.kT() = hf.kT();
    f.n() = hf.n();
    f.mu() = hf.mu();
    f.Vext() = 0;
    double energy = f.get_volume()*hf.energy();
    retval += check_functional_value("SFMTFluid", f, energy, 2e-11);
    //retval += f.run_finite_difference_test("SFMTFluid");
    Vector z = f.get_rz();
    for (int i=0; i<z.get_size(); i++) {
      if (z[i] > 0) {
        f.n()[i] *= 0.0001;
      }
    }
    retval += f.run_finite_difference_test("SFMTFluid inhomogeneous");
  }

  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
