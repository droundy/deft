#include "Monte-Carlo/square-well.h"

int num_errors = 0;

void check_dos(sw_simulation *sw) {
  double *D = new double[sw->energy_levels]();
  double *TD = new double[sw->energy_levels]();
  for (int i=0;i<sw->energy_levels;i++) {
    D[i] = exp(-sw->ln_energy_weights[i]);
  }

  for (int i=0;i<sw->energy_levels;i++) {
    double norm = 0;
    for (int de=-sw->biggest_energy_transition; de<=sw->biggest_energy_transition; de++) {
      norm += sw->transitions(i, de);
    }
    for (int j=0;j<sw->energy_levels;j++) {
      if (abs(i-j) <= sw->biggest_energy_transition) {
        TD[j] += D[i]*sw->transitions(i, j-i)/norm;
      }
    }
  }
  for (int i=0;i<sw->energy_levels;i++) {
    printf("D[%d] = %10g TD[%d] = %10g\n", i, D[i], i, TD[i]);
    double error = (D[i] - TD[i])/(D[i] + TD[i]);
    if (fabs(error) > 1e-10) {
      printf("FAIL error of %g at i=%d\n", error, i);
      num_errors++;
    }
  }
  delete[] D;
  delete[] TD;
}

int main(int argc, char **argv) {
  sw_simulation sw;
  sw.fractional_dos_precision = 1e-12;
  sw.energy_levels = 10;
  sw.min_energy_state = sw.energy_levels-1;
  sw.max_entropy_state = 0;

  sw.energy_histogram = new long[sw.energy_levels]();
  for (int i=0;i<sw.energy_levels;i++) sw.energy_histogram[i] = 1;
  sw.ln_energy_weights = new double[sw.energy_levels]();

  // Transitions from one energy to another
  sw.biggest_energy_transition = 10;
  sw.transitions_table =
    new long[sw.energy_levels*(2*sw.biggest_energy_transition+1)]();

  printf("\n** Flat density of states **\n\n");
  sw.transitions(0, 1) = 5;
  sw.transitions(0, 0) = 10;
  sw.transitions(sw.energy_levels-1, -1) = 5;
  sw.transitions(sw.energy_levels-1, 0) = 10;
  for (int i=1;i<sw.energy_levels-1;i++) {
    long counts = (73*long(i) + 37) % 1001; // pseudorandom
    sw.transitions(i,  0) = counts;
    sw.transitions(i,  1) = counts;
    sw.transitions(i, -1) = counts;
  }
  sw.update_weights_using_transitions();
  check_dos(&sw);
  for (int i=2;i<sw.energy_levels-1;i++) {
    const double dos0 = exp(sw.ln_energy_weights[1]);
    double error = (exp(sw.ln_energy_weights[i]) - dos0) / dos0;
    if (fabs(error) > 1e-10) {
      printf("FAIL error of %g at i=%d\n", error, i);
      num_errors++;
    }
  }

  printf("\n** Exponential density of states **\n\n");
  for (int i=1;i<sw.energy_levels-1;i++) {
    long counts = (73*long(i) + 37) % 1001; // pseudorandom
    sw.transitions(i,  0) = counts;
    sw.transitions(i,  1) = i*counts;
    sw.transitions(i, -1) = counts;
  }
  sw.update_weights_using_transitions();
  check_dos(&sw);

  printf("\n** Beyond exponential density of states **\n\n");
  for (int i=1;i<sw.energy_levels-1;i++) {
    sw.transitions(i,  0) = 10;
    sw.transitions(i,  1) = 500;
    sw.transitions(i, -1) = 1;
  }
  sw.update_weights_using_transitions();
  check_dos(&sw);

  return num_errors;
}
