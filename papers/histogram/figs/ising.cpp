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
double iterations = 0; 
double J = 0;  // Define our spin coupling.
const int N = 25;  // Define number of spin sites.
// user should be able to define minT? 

// ---------------------------------------------------------------------
// Mersenne Twister Engine -- build the random seed generator
// ---------------------------------------------------------------------

//the random device that will seed the generator
std::random_device seeder;
// The Mersenne twister engine which has a distinct advantage
// over rand().  Only for c++11 compilers. 
std::mt19937 engine(seeder());
std::uniform_int_distribution<int> spin_dist(-1, 1);

// ---------------------------------------------------------------------
// Generate 2D Spin Lattice 
// ---------------------------------------------------------------------

int spin_lattice[N][N];
int spin_number = spin_dist(engine); // currently either -1,0,+1

//inline generate_spin_lattice() {
  //for (int i = 0; i < 2; i++) {
    //for (int j = 0; j < 2; j++) {
      //spin_lattice[i][j] = spin_number;
    //}
  //}
  //return spin_lattice;
//}

// ---------------------------------------------------------------------
// Flip a spin site randomly
// ---------------------------------------------------------------------

std::uniform_int_distribution<int> X(0, N);
std::uniform_int_distribution<int> Y(0, N);

//int spin_site_flip[X][Y];

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
  printf("lets randomize: %d\n", spin_number);
}

