#include "Minimizer.h"
#include <stdio.h>

class SqrSum : public Functional {
  double energy(const Vector &x, Verbosity v) const {
    const int sz = x.get_size();
    double out = 0;
    for (int i=0; i<sz; i++) {
      if (v >= verbose) printf("x[%2d] = %g\n", i, x[i]);
      out += double(i+0.1)*x[i]*x[i];
    }
    return out;
  }
  Vector grad(const Vector &x, const Bools *) const {
    return 2*x;
  }
};

int main() {
  int errors = 0;
  const int N = 10;

  {
    printf("*** Testing that minimization gets (almost) exact answer in easy case ***\n");
    Vector foo(N);
    for (int i=0;i<N;i++) {
      foo[i] = i*0.1;
    }
    SqrSum sqr;
    Minimizer min(&sqr, &foo);
    printf("Starting energy is %g\n\n", min.energy());
    while (min.improve_energy(chatty)) {
    }
    min.print_info();
    // We should be able to quickly find the exact minimum, which
    // happens to be zero...
    if (min.energy() > 1e-320) {
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
    while (min.improve_energy(chatty)) {
    }
    min.print_info();
    // We should be able to quickly find the exact minimum, which
    // happens to be zero...
    if (min.energy() > prec) {
      printf("FAIL: Energy error is too big! %g vs %g\n", min.energy(), prec);
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
