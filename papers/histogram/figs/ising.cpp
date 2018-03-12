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

double E = 0;  // Define our initial system energy.
double M = 0;  // Define our initial system magnetization.
double J = 0;  // Define our spin coupling.

const int N = 25;  // Define number of spin sites.
const int Q = 2;   // Define number of spin states.

double minT = 0; // Define minimum Temperature.

// ---------------------------------------------------------------------
// ising.h files -- this could likely go into an ising.h file
// ---------------------------------------------------------------------

struct ising_simulation {
  long iteration;   // the current iteration number

  int spin_number;        // the current spin number
  int site_numberI;
  int site_numberJ;
  int spin_lattice[N][N]; // the size of the spin system
  
  void flip_a_spin(); // attempt to flip a spin
  void initialize_spin_lattice(); // generate the spin lattice
  void random_site_flip(); // flip the spin of a random site
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
    spin_number = spin_number - (Q + 1);
  }
}

// ---------------------------------------------------------------------
// Generate 2D Spin Lattice 
// ---------------------------------------------------------------------

void ising_simulation::initialize_spin_lattice() {
  assert (Q%2 == 0);  // terminate program if Q is odd!

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      flip_a_spin();
      spin_lattice[i][j] = spin_number;
    }
  }
}

// ---------------------------------------------------------------------
// Flip a spin site randomly
// ---------------------------------------------------------------------

void ising_simulation::random_site_flip() {
    std::random_device seeder;
    std::mt19937 engine(seeder());
    std::uniform_int_distribution<int> site_distI(1, N);
    std::uniform_int_distribution<int> site_distJ(1, N);
    
    site_numberI = site_distI(engine); // currently ranges from 1 --> N
    site_numberJ = site_distJ(engine); // currently ranges from 1 --> N
    
    flip_a_spin();

    spin_lattice[site_numberI][site_numberJ] = spin_number;
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
  took("Starting program");
  printf("version: %s\n",version_identifier());
}

