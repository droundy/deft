#include <stdio.h>
#include <float.h>
#include <time.h>
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <popt.h>
#include <sys/stat.h>
#include <string.h>
#include <errno.h>
#include <popt.h>
#include "handymath.h"
#include "vector3d.h"

#include "version-identifier.h"

// ---------------------------------------------------------------------
// Define "Constants" -- set from arguments then unchanged
// ---------------------------------------------------------------------

double M = 0;      // Define our initial system magnetization.
double J = 1;      // Define our spin coupling.
double minT = 0.2; // Define minimum Temperature.

const int Q = 2;   // Define number of spin states.

// ---------------------------------------------------------------------
// ising.h files -- this could likely go into an ising.h file
// ---------------------------------------------------------------------

enum dos_types { histogram_dos, transition_dos, weights_dos };

struct ising_simulation {
  long iteration;   // the current iteration number
  int N;            // N*N is the number of sites
  int *S;
  double E;      // system energy.
  double T;      // canonical temperature.

  long energy_levels;
  long energies_found;
  long *energy_histogram;

  explicit ising_simulation(int N); // generate the spin lattice

  int random_flip(int oldspin) const; // pick a new random (but changed) spin

  void reset_histograms();
  void flip_a_spin();
  double calculate_energy();
  
  double* compute_ln_dos(dos_types dos_type);
  double *ln_energy_weights;
  double *ln_dos;
  int max_entropy_state;
};

// ising_simulation methods

void ising_simulation::reset_histograms(){

  for(int i = 0; i < energy_levels; i++){
    energy_histogram[i] = 0;
  }
}

int ising_simulation::random_flip(int oldspin) const {
  if (Q==2) return -oldspin;

  return random::ran64() % Q - Q/2; // gives value -Q/2 to Q/2-1?
}

ising_simulation::ising_simulation(int NN) {
  iteration = 0;
  E = 0;
  T = 0;
  N = NN;
  S = new int[N*N];
  // Energy histogram
  // abs. ground state energy for 2D Ising model.  But with new long rather
  // than just long for energy_histogram could I keep track of energy levels
  // found and update as needed?
  // i.e. if H(E) == 0 then energy_levels ++ "we have found a new energy".
  energy_levels = 2*J*N*N;
  ln_energy_weights = new double[energy_levels]();
  energy_histogram = new long[energy_levels]();
  max_entropy_state = 0;
  //printf("energy_levels %ld\n", energy_levels);
  reset_histograms();

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      S[i+N*j] = 1; // initialize to all 1s
    }
  }
}

void ising_simulation::flip_a_spin() {
  iteration += 1;
  int k = random::ran64() % (N*N);
  int i = k % N;
  int j = k / N;
  int old = S[j + i*N];
  S[j + i*N] = random_flip(S[j + i*N]);

  int neighbor_spins = S[j + ((i+1) % N)*N] + S[j + ((i-1+N) % N)*N] +
                       S[((j+1) % N) + i*N] + S[((j-1+N) % N) + i*N];

// ---------------------------------------------------------------------
// Canonical Monte-Carlo
// ---------------------------------------------------------------------
  const double deltaE = J*(old - S[j + i*N])*neighbor_spins;
  const double lnprob = -deltaE/T;

  if (lnprob < 0 && random::ran() > exp(lnprob)) {
    S[j + i*N] = old;
    //printf("not flipping from E=%g\n", E);
  } else {
    E += deltaE;

//  if (energy_histogram[abs(lround(E))] == 0) {
//    energies_found++;  // we have found a new energy!
//  };
    //printf("flipping gives E=%g\n", E);
  }

  energy_histogram[abs(lround(E))] += 1;
  //energy_histogram[energies_found] += 1;

// ---------------------------------------------------------------------
// Broad-Histogram Methods
// ---------------------------------------------------------------------
}

double ising_simulation::calculate_energy() {
  E = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // we want to enforce periodic boundary conditions.
      int neighbor_spins = S[(i+1) % N + j*N] + S[i + ((j+1) % N)*N];

      E += -J*neighbor_spins*S[i+j*N];
    }
  }
  return E;
}

double* ising_simulation::compute_ln_dos(dos_types dos_type) {

  double *ln_dos = new double[energy_levels]();

  if(dos_type == histogram_dos){
    for(int i = max_entropy_state; i < energy_levels; i++){
      if(energy_histogram[i] != 0){
        ln_dos[i] = log(energy_histogram[i]) - ln_energy_weights[i];
      }
      else ln_dos[i] = -DBL_MAX; //apparently in float.h?
    }
  }
  return ln_dos;
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

  poptContext optCon;

  // -------------------------------------------------------------------
  // Parse input options
  // -------------------------------------------------------------------

  poptOption optionsTable[] = {

    /*** ISING MODEL PARAMETERS ***/

    //{"N", '\0', POPT_ARG_INT, &ising.N, 0, "N*N is the number of spin sites", "INT"},
    //{"Q", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &ising.Q, 0,
     //"The number of spin states", "INT"},

    /*** SIMULATION ITERATIONS ***/

    //{"iterations", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &simulation_iterations,
     //0, "Number of iterations for which to run the simulation", "INT"},
    //{"round_trips", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &simulation_round_trips,
     //0, "Number of round trips (pessimistic samples) to run the simulation", "INT"},

    /*** PARAMETERS DETERMINING OUTPUT FILE DIRECTORY AND NAMES ***/

    //{"data_dir", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &data_dir, 0,
     //"Directory in which to save data", "data_dir"},
    //{"filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     //"Base of output file names", "STRING"},
    //{"filename_suffix", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
     //&filename_suffix, 0, "Output file name suffix", "STRING"},

    /*** HISTOGRAM METHOD OPTIONS ***/

    //{"kT", '\0', POPT_ARG_DOUBLE, &fix_kT, 0, "Use a fixed temperature of kT"
     //" rather than adjusted weights", "DOUBLE"},

    /*** HISTOGRAM METHOD PARAMETERS ***/

    //{"wl_factor", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &wl_factor,
     //0, "Initial value of Wang-Landau factor", "DOUBLE"},

    /*** END CONDITION PARAMETERS ***/

    //{"min_samples", '\0', POPT_ARG_INT, &sw.min_samples, 0,
     //"Number of times to sample mininum energy", "INT"},
    //{"init_iters", '\0', POPT_ARG_INT, &sw.init_iters, 0,
     //"Number of iterations for initialization", "INT"},
    //{"min_T", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     //&ising.min_T, 0, "The minimum temperature that we care about", "DOUBLE"},

    POPT_AUTOHELP
    POPT_TABLEEND
  };
  optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTION...]\nRequired arguments: square-root of number of sites (N), "
                         "and the number of spin states (Q)");

  int c = 0;
  // go through arguments, set them based on optionsTable
  while((c = poptGetNextOpt(optCon)) >= 0);
  if (c < -1) {
    fprintf(stderr, "\n%s: %s\n", poptBadOption(optCon, 0), poptStrerror(c));
    return 1;
  }
  poptFreeContext(optCon);

  printf("------------------------------------------------------------------\n");
  printf("Running %s with parameters:\n", argv[0]);
  for(int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') printf("\n");
    printf("%s ", argv[i]);
  }

  ising_simulation ising(4);
  ising.T = 2;
  ising.energies_found = 0; // we haven't found any energies yet.

  took("Starting program");
  printf("version: %s\n",version_identifier());

  long total_iteration = 12000;

  ising.calculate_energy();
  for (long i = 0; i < total_iteration; i++) {
    ising.flip_a_spin();
  }
  printf("I think the energy is %g\n", ising.E);
  ising.calculate_energy();
  printf("the energy is %g\n", ising.E);

  took("Running");

  ising.ln_dos = ising.compute_ln_dos(histogram_dos);

  for(int i = 0; i <= ising.energy_levels; i++){
    printf("histogram is %ld\n while lndos is %g\n", ising.energy_histogram[i], ising.ln_dos[i]);
  }
  //printf("lndos is %f\n", ising.ln_dos[1]);
  for(int i = 0; i <= ising.energy_levels; i++){
    if (ising.energy_histogram[i] > 0) {
      ising.energies_found++;
      printf("energy histogram at %ld is %ld\n",
              ising.energies_found,ising.energy_histogram[i]);
    }
  }
  // -------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // -------------------------------------------------------------------

  delete[] ising.energy_histogram;

}

