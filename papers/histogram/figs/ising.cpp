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
#include <iostream>
#include "version-identifier.h"

// ---------------------------------------------------------------------
// Define "Constants" -- set from arguments then unchanged
// ---------------------------------------------------------------------

double M = 0;      // Define our initial system magnetization.
int J = 1;      // Define our spin coupling.
double minT = 0.2; // Define minimum Temperature.

const int Q = 2;   // Define number of spin states.

// ---------------------------------------------------------------------
// ising.h files -- this could likely go into an ising.h file
// ---------------------------------------------------------------------

enum dos_types { histogram_dos, transition_dos, weights_dos };

struct ising_simulation {
  long moves;   // the current number of moves
  int N;            // N*N is the number of sites
  int *S;
  int E;      // system energy.

  void canonical_temperature(double t);
  double T;      // canonical temperature.

  long energy_levels;
  long index_from_energy(long E) const {
    return J*N*N/2 + E/4;
  }

  long energies_found;
  long *energy_histogram;

  explicit ising_simulation(int N); // generate the spin lattice
  ~ising_simulation();

  int random_flip(int oldspin) const; // pick a new random (but changed) spin

  void flip_a_spin();
  double calculate_energy();

  double* compute_ln_dos(dos_types dos_type);
  double *ln_energy_weights;
  double *ln_dos;
  int max_entropy_state;
};

// ising_simulation methods

int ising_simulation::random_flip(int oldspin) const {
  if (Q==2) return -oldspin;

  return random::ran64() % Q - Q/2; // gives value -Q/2 to Q/2-1?
}

ising_simulation::ising_simulation(int NN) {
  moves = 0;
  E = 0;
  T = 0;
  N = NN;
  S = new int[N*N];
  // energy histogram
  energy_levels = J*N*N;
  ln_energy_weights = new double[energy_levels]();
  energy_histogram = new long[energy_levels]();
  max_entropy_state = 0;
  energies_found = 0; // we haven't found any energies yet.
  //printf("energy_levels %ld\n", energy_levels);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      S[i+N*j] = 1; // initialize to all 1s
    }
  }
}

ising_simulation::~ising_simulation() {
  delete[] energy_histogram;
  delete[] ln_energy_weights;
  delete[] S;
}

void ising_simulation::canonical_temperature(double t) {
  T = t;
  for (int i=0; i<energy_levels; i++) {
    ln_energy_weights[i] = 4*i/T;
  }
}

void ising_simulation::flip_a_spin() {
  moves += 1;
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
  const int deltaE = J*(old - S[j + i*N])*neighbor_spins;
  const double lnprob = -deltaE/T; // ln_energy_weights[-E] - ln_energy_weights[-(E+deltaE)];

  //printf("compare %g with %g %d->%d\n", -deltaE/T,
        //ln_energy_weights[index_from_energy(E)]
         //- ln_energy_weights[index_from_energy(E+deltaE)],
        //E, E+deltaE);

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

  energy_histogram[index_from_energy(E)] += 1;

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

  if (dos_type == histogram_dos) {
    for (int i = max_entropy_state; i < energy_levels; i++) {
      if (energy_histogram[i] != 0) {
        ln_dos[i] = log(energy_histogram[i]) + ln_energy_weights[i];
      } else {
        ln_dos[i] = -DBL_MAX; // located in <float.h>.
      }
    }
  }
  return ln_dos;
}

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

// ---------------------------------------------------------------------
// Initialize Main
// ---------------------------------------------------------------------

int main(int argc, const char *argv[]) {
  random::seed(0);
  // some miscellaneous default or dummy simulation parameters

  bool Jordan = true; 
  int NN = 10;
  int resume = false;
  long total_moves = 10000;

  char *data_dir = new char[1024];
  sprintf(data_dir,"none");
  char *default_data_dir = new char[1024];
  sprintf(default_data_dir, "papers/histogram/data/ising");
  char *filename = new char[1024];
  sprintf(filename, "none");
  char *filename_suffix = new char[1024];
  sprintf(filename_suffix, "none");

  //int Q = 2;  // defualt to ising model
  poptContext optCon;

  // -------------------------------------------------------------------
  // Parse input options
  // -------------------------------------------------------------------

  poptOption optionsTable[] = {
    {"resume", '\0', POPT_ARG_NONE, &resume, 0,
     "Resume previous simulation", "BOOLEAN"},
    //{"seed", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0,
    // "Seed for the random number generator", "INT"},

    /*** ISING MODEL PARAMETERS ***/

    {"N", '\0', POPT_ARG_INT, &NN, 0, "N*N is the number of spin sites", "INT"},
    //{"Q", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &ising.Q, 0,
    // "The number of spin states", "INT"},

    /*** SIMULATION ITERATIONS ***/

    {"total-moves", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &total_moves,
     0, "Number of moves for which to run the simulation", "INT"},

    /*** PARAMETERS DETERMINING OUTPUT FILE DIRECTORY AND NAMES ***/

    {"dir", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &data_dir, 0,
     "Directory in which to save data", "data_dir"},
    {"filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     "Base of output file names", "STRING"},
     

    /*** HISTOGRAM METHOD OPTIONS ***/

    //{"kT", '\0', POPT_ARG_DOUBLE, &fix_kT, 0, "Use a fixed temperature of kT"
    // " rather than adjusted weights", "DOUBLE"},

    /*** HISTOGRAM METHOD PARAMETERS ***/

    //{"wl_factor", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &wl_factor,
     //0, "Initial value of Wang-Landau factor", "DOUBLE"},

    /*** END CONDITION PARAMETERS ***/

    //{"min_samples", '\0', POPT_ARG_INT, &sw.min_samples, 0,
     //"Number of times to sample mininum energy", "INT"},

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

  ising_simulation ising(NN);
  ising.canonical_temperature(6);

  took("Starting program");
  printf("version: %s\n",version_identifier());

  ising.calculate_energy();
  for (long i = 0; i < total_moves; i++) {
    ising.flip_a_spin();
  }

  printf("I think the energy is %d\n", ising.E);
  ising.calculate_energy();
  printf("the energy is %d\n", ising.E);

  took("Running");

  ising.ln_dos = ising.compute_ln_dos(histogram_dos);

//  for(int i = 0; i <= ising.energy_levels; i++){
//    printf("histogram is %ld\n while lndos is %g\n", ising.energy_histogram[i], ising.ln_dos[i]);
//  }

  for(int i = 0; i <= ising.energy_levels; i++){
    if (ising.energy_histogram[i] > 0) {
      ising.energies_found++;
      printf("energy histogram at %ld is %ld with lndos %g\n",
              ising.energies_found,ising.energy_histogram[i],ising.ln_dos[i]);
    }
  }

  // ----------------------------------------------------------------------------
  // Generate save file info
  // ----------------------------------------------------------------------------

  // Set default data directory
  if (strcmp(data_dir,"none") == 0) {
    sprintf(data_dir,"%s",default_data_dir);
    printf("\nUsing default data directory: [deft]/%s\n",data_dir);
  }

  mkdir(data_dir, 0777); // create save directory

  char *dos_fname = new char[1024];
  sprintf(dos_fname, "%s/%s-dos.dat", data_dir, filename);

  //if (fix_kT) {
  //  sprintf(headerinfo,
  //          "%s# histogram method: canonical (fixed temperature)\n"
  //          "# kT: %g\n",
  //          headerinfo, fix_kT);
  //}

  // ----------------------------------------------------------------------------
  // Resume functionality
  // ----------------------------------------------------------------------------

  if (resume) {
     // We are continuing a previous simulation.
    if (Jordan) {// This will be a Monte Carlo method eventually.
      // Open the file!
      FILE *rfile = fopen((const char *)dos_fname, "r");
      if (rfile != NULL) {
        printf("I'm resuming with the method.\n");
        fscanf(rfile, "import numpy as np\n");
        random::resume_from_dump(rfile);
        char * line = new char[1000];
        if (fscanf(rfile, " %[^\n]\n", line) != 1) {
          printf("error reading headerinfo!\n");
          exit(1);
        }
        printf("line: '%s'\n",line);
        if (fscanf(rfile, " NN = %d\n", &NN) != 1) {
          printf("error reading NN!\n");
          exit(1);
        }
        printf("NN is now %d\n", NN);
        if (fscanf(rfile, " total-moves = %ld\n", &total_moves) != 1) {
          printf("error reading total-moves!\n");
          exit(1);
        }
        printf("total-moves is now %ld\n", total_moves);
        fclose(rfile);
      }
      // Need to read and set: iterations, ln_energy_weights (from ln_dos)
      } else {
        printf("I do not know how to resume yet!\n");
        exit(1);
      }
  }


  // Save energy histogram
      {
        FILE *dos_out = fopen((const char *)dos_fname, "w");

        fprintf(dos_out, "import numpy as np\n\n");
        random::dump_resume_info(dos_out);
        fprintf(dos_out,"version = %s\n\n",version_identifier());

        fprintf(dos_out,"NN = %i\n",NN);
        fprintf(dos_out,"total-moves = %li\n",total_moves);
        fprintf(dos_out,"minT = %g\n",minT);
        fprintf(dos_out,"moves = %li\n",ising.moves);
        fprintf(dos_out,"E = %i\n",ising.E);
        fprintf(dos_out,"J = %i\n", J);
        // inserting arrays into text file.
        fprintf(dos_out, "lndos = np.array([\n");
        for (int i = 0; i < ising.energy_levels; i++) {
          fprintf(dos_out,"\t%.17g,\n", ising.ln_dos[i]);
        }
        fprintf(dos_out, "])\n");
        fprintf(dos_out, "lnw = np.array([\n");
        for (int i = 0; i < ising.energy_levels; i++) {
          fprintf(dos_out,"\t%.17g,\n", ising.ln_energy_weights[i]);
        }
        fprintf(dos_out, "])\n");
        fprintf(dos_out, "histogram = np.array([\n");
        for (int i = 0; i < ising.energy_levels; i++) {
          fprintf(dos_out,"\t%ld,\n", ising.energy_histogram[i]);
        }
        fprintf(dos_out, "])\n");
        fprintf(dos_out, "S = np.array([\n");
        for (int i = 0; i < NN*NN; i++) {
          fprintf(dos_out,"\t%d,\n", ising.S[i]);
        }
        fprintf(dos_out, "])\n");

        fclose(dos_out);
      }
  // -------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // -------------------------------------------------------------------
}

