#include <stdio.h>
#include <time.h>
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <popt.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include "handymath.h"
#include "vector3d.h"
#include "random"

#include "version-identifier.h"

// ---------------------------------------------------------------------
// Define "Constants" -- set from arguments then unchanged
// ---------------------------------------------------------------------

double E = 0;      // Define our initial system energy.
double M = 0;      // Define our initial system magnetization.
double Jc = 1;      // Define our spin coupling.
double minT = 0.2; // Define minimum Temperature.

const int N = 20;  // Define number of spin sites.
const int Q = 2;   // Define number of spin states.

// ---------------------------------------------------------------------
// ising.h files -- this could likely go into an ising.h file
// ---------------------------------------------------------------------

struct ising_simulation {
  long iteration;   // the current iteration number

  int spin_number;        // the current spin number
  int spin_number_at_site[N*N];
  int I;                  // random row coord spin number
  int J;                  // random col coord spin number
  int spin_lattice[N][N]; // the size of the spin system
  int neighbor_spins;

  double prob;

  void flip_a_spin(); // attempt to flip a spin
  int* initialize_spin_lattice(); // generate the spin lattice

  int* canonical_flip();
  double calculate_energy();
};

// ---------------------------------------------------------------------
// Mersenne Twister Engine -- build the random seed generator
// ---------------------------------------------------------------------

// ising_simulation methods

void ising_simulation::flip_a_spin() {
  //the random device that will seed the generator
  std::random_device seeder;
  // The Mersenne twister engine which has a distinct advantage
  // over rand().  Only for c++11 compilers.
  std::mt19937 engine(seeder());
  std::uniform_int_distribution<int> spin_dist(1, Q);

  spin_number = spin_dist(engine); // currently ranges from 1 --> Q

  if (spin_number > Q/2.0) {
    spin_number += - Q - 1;
  }
}

// ---------------------------------------------------------------------
// Generate 2D Spin Lattice
// ---------------------------------------------------------------------

// I want to make this a pointer so I can calculate energy from each
// spin number at a lattice site.
int* ising_simulation::initialize_spin_lattice() {
  assert (Q%2 == 0);  // terminate program if Q is odd!

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      flip_a_spin();
      spin_lattice[i][j] = spin_number;
      spin_number_at_site[i*j] = spin_lattice[i][j];
    }
  } return spin_number_at_site;
}

// ---------------------------------------------------------------------
// Canonical Monte-Carlo
// ---------------------------------------------------------------------

int* ising_simulation::canonical_flip() {

  int i;
  int j;

  std::random_device seeder;
  std::mt19937 engine(seeder());
  std::uniform_int_distribution<int> site_distI(1, N);
  std::uniform_int_distribution<int> site_distJ(1, N);

  I = site_distI(engine); // currently ranges from 1 --> N
  J = site_distJ(engine); // currently ranges from 1 --> N

  flip_a_spin();
  //printf("spin number is really %d\n", spin_number);
  spin_lattice[i][j] = spin_number;

  // we want to enforce periodic boundary conditions.
  i = (I + N) % N; // coerce in range: takes care of wrapping
  j = (J + N) % N;
  neighbor_spins = spin_lattice[i+1][j] + spin_lattice[i-1][j] +
                   spin_lattice[i][j+1] + spin_lattice[i][j-1];

  prob = 2*spin_lattice[I][J]*neighbor_spins;

  if (prob < 0) {
    spin_lattice[I][J] *= -1;
  } else if (rand() < exp(-prob*minT)) {
    spin_lattice[I][J] *= -1;
  }
  return spin_number_at_site;
}

double ising_simulation::calculate_energy() {

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // we want to enforce periodic boundary conditions.
      I = (i + N) % N; // coerce in range: takes care of wrapping
      J = (j + N) % N;
      neighbor_spins = spin_lattice[I+1][J] + spin_lattice[I-1][J] +
                       spin_lattice[I][J+1] + spin_lattice[I][J-1];
      
      printf("spin number is %d\n at (%d,%d)\n", spin_lattice[I+1][J], I+1,J);
      E += -Jc*neighbor_spins*spin_lattice[i][j];
    }
  } return E;
}
// ---------------------------------------------------------------------
// Initialize Main
// ---------------------------------------------------------------------

static double took(const char *name) {
  assert(name); // so it'll count as being used...
  static clock_t last_time = clock();
  clock_t t = clock();
  double seconds = (t-last_time)/double(CLOCKS_PER_SEC);
  if (seconds > 120) {
    printf("%s took %.0f minutes and %g seconds.\n", name, seconds/60,
           fmod(seconds,60));
  } else {
    printf("%s took %g seconds...\n", name, seconds);
  }
  fflush(stdout);
  last_time = t;
  // We round our times to a close power of e to make our code more
  // likely to be reproducible, i.e. so two runs with identical
  // parameters (including random number seed) should (unless we are
  // unlucky) result in identical output.
  return exp(ceil(log(seconds)));
}

int main(int argc, const char *argv[]) {

  ising_simulation ising;

  took("Starting program");
  printf("version: %s\n",version_identifier());

  long iteration = 500;
  ising.initialize_spin_lattice();

  for (int i = 0; i < iteration; i++) {
    ising.canonical_flip();
  }
  ising.calculate_energy();
  printf("the energy is %g\n", E/iteration);

}

