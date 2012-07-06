#include "Minimizer.h"
#include <stdio.h>

inline double spring(int i) {
  return double(i+1);
  return 1.0;
}

class SqrSum : public Functional {
  double energy(const Vector &x) const {
    const int sz = x.get_size();
    double out = 0;
    for (int i=0; i<sz; i++) {
      out += spring(i)*x[i]*x[i];
    }
    return out;
  }
  double energy_per_volume(const Vector &x) const {
    return 0;
  }
  double denergy_per_volume_dx(const Vector &x) const {
    return 0;
  }
  Vector grad(const Vector &x) const {
    const int sz = x.get_size();
    Vector out(sz);
    for (int i=0; i<sz; i++) {
      out[i] = spring(i)*2*x[i];
    }
    return out;
  }
  EnergyGradAndPrecond energy_grad_and_precond(const Vector &x) const {
    EnergyGradAndPrecond egpg;
    egpg.energy = energy(x);
    egpg.grad = grad(x);
    egpg.precond = grad(x);
    const int sz = x.get_size();
    for (int i=0; i<sz; i++) {
      egpg.precond[i] /= spring(i);
    }
    return egpg;
  }
  void printme(const char *) const {
  }
  bool have_preconditioner() const { return true; }
};

int main() {
  int errors = 0;
  const int N = 1000;

  {
    printf("*** Testing that minimization gets (almost) exact answer in easy case ***\n");
    Vector foo(N);
    for (int i=0;i<N;i++) {
      foo[i] = i*0.1;
    }
    SqrSum sqr;
    Minimizer min(&sqr, &foo);
    printf("Starting energy is %g\n\n", min.energy());
    while (min.improve_energy(quiet)) {
    }
    min.print_info();
    // We should be able to quickly find the exact minimum, which
    // happens to be zero...
    if (min.energy() > 1e-300) {
      printf("FAIL: Energy is too big! %g\n", min.energy());
      errors++;
    }
  }

  {
    printf("*** Testing that minimization can get approximate answer in easy case ***\n");
    Vector foo(N);
    for (int i=0;i<N;i++) {
      foo[i] = i*0.1;
    }
    SqrSum sqr;
    Minimizer min(&sqr, &foo);
    const double prec = 1e-9;
    min.set_precision(prec);
    printf("Starting energy is %g\n\n", min.energy());
    while (min.improve_energy(quiet)) {
    }
    min.print_info();
    // We should be able to quickly find the exact minimum, which
    // happens to be zero...
    if (min.energy() > prec) {
      printf("FAIL: Energy error is too big! %g vs %g\n", min.energy(), prec);
      errors++;
    }
    const double condition_number = N; // I'm not sure about this...
    printf("Took %d iterations, and I expected no more than %g\n",
           min.get_iteration_count(), condition_number);
    if (min.get_iteration_count() > condition_number) {
      printf("FAIL: Took too many iterations! %d vs %g\n", min.get_iteration_count(), condition_number);
      errors++;
    }
  }

  {
    printf("\n*** Testing minimization with preconditioning ***\n");
    Vector foo(N);
    for (int i=0;i<N;i++) {
      foo[i] = i*0.1;
    }
    SqrSum sqr;
    Minimizer min(&sqr, &foo);
    const double prec = 1e-9;
    min.set_precision(prec);
    min.precondition(true);
    printf("Starting energy is %g\n\n", min.energy());
    while (min.improve_energy(louder(min_details))) {
    }
    min.print_info();
    // We should be able to quickly find the exact minimum, which
    // happens to be zero...
    if (min.energy() > prec) {
      printf("FAIL: Energy error is too big! %g vs %g\n", min.energy(), prec);
      errors++;
    }
    const double condition_number = 2; // Preconditioning is exact in this case
    printf("Took %d iterations, and I expected no more than %g\n",
           min.get_iteration_count(), condition_number);
    if (min.get_iteration_count() > condition_number) {
      printf("FAIL: Took too many iterations! %d vs %g\n", min.get_iteration_count(), condition_number);
      errors++;
    }
  }

  if (errors == 1) {
    printf("There was %d error!\n", errors);
  } else if (errors > 1) {
    printf("There were %d errors!\n", errors);
  } else {
    printf("All good!\n");
  }
  return errors;
}
