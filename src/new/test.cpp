#include "Minimizer.h"
#include <stdio.h>

class SqrSum : public Functional {
  double energy(const Vector &x, Verbosity v) const {
    const int sz = x.get_size();
    double out = 0;
    for (int i=0; i<sz; i++) {
      if (v >= verbose) printf("x[%2d] = %g\n", i, x[i]);
      out += x[i]*x[i];
    }
    return out;
  }
  Vector grad(const Vector &x, const Bools *) const {
    return 2*x;
  }
};

int main() {
  const int N = 10;
  Vector foo(N);
  for (int i=0;i<N;i++) {
    foo[i] = i*0.1;
  }
  SqrSum sqr;
  Minimizer min(&sqr, &foo);
  min.set_precision(1e-9);
  printf("Starting energy is %g\n\n", min.energy());
  while (min.improve_energy(chatty)) {
  }
  min.print_info();
  // We should be able to quickly find the exact minimum, which
  // happens to be zero...
  if (min.energy() != 0.0) {
    printf("Energy is too big! %g\n", min.energy());
    return 1;
  }
  printf("All good!\n");
  return 0;
}
