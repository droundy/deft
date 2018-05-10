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

  /* The following determine how we perform an ising flip */
  double wl_factor; // if it is non-zero, update weights on each move WL-style.
  double sa_t0; // if it is non-zero, update wl_factor on each move SA-style.
  double sa_prefactor; // prefactor in computing wl_factor when running SA.
  int use_sad;  // if nonzero, dynamically update sa_t0 when we
                // encounter new energies, value is the fraction to
                // use.
  bool use_wl; // if true, we are using WL

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
  void end_flip_updates();
  double calculate_energy();

  void compute_ln_dos(dos_types dos_type);
  double *ln_energy_weights;
  double *ln_dos;
  int max_entropy_state, min_energy_state, min_important_energy;
  int too_high_energy; // too_high_energy is used in SAD, and is the
                       // highest ever value that max_entropy_state has
                       // taken.  It is used to define the range over
                       // which our histogram is made flat.
  int too_low_energy; // too_low_energy is used in SAD, and is the
                      // lowest ever value that min_important_energy has
                      // taken.  It is used to define the range over
                      // which our histogram is made flat.

  /*** HISTOGRAM METHODS ***/

};

// ising_simulation methods

int ising_simulation::random_flip(int oldspin) const {
  if (Q==2) return -oldspin;

  return random::ran64() % Q - Q/2; // gives value -Q/2 to Q/2-1?
}

ising_simulation::ising_simulation(int NN) {
  use_sad = 0; // default to not using SAD MC (think happy thoughts!)
  moves = 0;
  E = 0;
  T = 0;
  N = NN;
  S = new int[N*N];
  // energy histogram
  energy_levels = J*N*N;
  ln_energy_weights = new double[energy_levels]();
  ln_dos = new double[energy_levels]();
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
  delete[] ln_dos;
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

  const int deltaE = J*(old - S[j + i*N])*neighbor_spins;

  double prob = 1;
  if (use_sad) {
    const int e1 = E;
    const int e2 = E + deltaE;
    double lnw1 = ln_energy_weights[e1];
    if (e1 < too_high_energy) {
      lnw1 = ln_energy_weights[too_high_energy];
    } else if (e1 > too_low_energy) {
      lnw1 = ln_energy_weights[too_low_energy] + (e1 - too_low_energy)/minT;
    }
    double lnw2 = ln_energy_weights[e2];
    if (e2 < too_high_energy) {
      lnw2 = ln_energy_weights[too_high_energy];
    } else if (e2 > too_low_energy) {
      lnw2 = ln_energy_weights[too_low_energy] + (e2 - too_low_energy)/minT;
    }
    const double lnprob = lnw2 - lnw1;
    prob = exp(lnprob); //WHAT TO DO HERE?
    //if (e2 < too_high_energy && e1 != e2) {
      //printf("prob hi %d -> %d = %g\n", e1, e2, Pmove);
    //} else if (e2 > too_low_energy && e1 != e2) {
      //printf("prob lo %d -> %d = %g\n", e1, e2, Pmove);
    //}
  } else {
    const double lnprob = ln_energy_weights[-E] - ln_energy_weights[-(E+deltaE)]; //-deltaE/T;
      if (lnprob < 0 && random::ran() > exp(lnprob)) {
        S[j + i*N] = old;
        //printf("not flipping from E=%g\n", E);
        } else {
          E += deltaE;
        }
    }
    if (prob < 1) {
      if (random::ran() > prob) {
        // We want to reject this move because it is too improbable
        // based on our weights.
        //if (get_new_neighbors) delete[] temp.neighbors;
          //end_move_updates();
          //return;
      }
    }

// ---------------------------------------------------------------------
// Canonical Monte-Carlo
// ---------------------------------------------------------------------
  //const int deltaE = J*(old - S[j + i*N])*neighbor_spins;
  //const double lnprob = -deltaE/T; // ln_energy_weights[-E] - ln_energy_weights[-(E+deltaE)];

  //printf("compare %g with %g %d->%d\n", -deltaE/T,
        //ln_energy_weights[index_from_energy(E)]
         //- ln_energy_weights[index_from_energy(E+deltaE)],
        //E, E+deltaE);

  //if (lnprob < 0 && random::ran() > exp(lnprob)) {
  //  S[j + i*N] = old;
  //  //printf("not flipping from E=%g\n", E);
  //} else {
  //  E += deltaE;

//  if (energy_histogram[abs(lround(E))] == 0) {
//    energies_found++;  // we have found a new energy!
//  };
    //printf("flipping gives E=%g\n", E);
  //}

  energy_histogram[index_from_energy(E)] += 1;
  // end move updates here???

// ---------------------------------------------------------------------
// Broad-Histogram Methods
// ---------------------------------------------------------------------
}

void ising_simulation::end_flip_updates(){
   // update iteration counter, energy histogram, and walker counters
  //if(moves % N*N == 0) iteration++; ???
  static int max_energy_seen = -1;
  static int min_energy_seen = -1;
  if (sa_t0 || use_sad) {
    if (energy_histogram[E] == 0) {
      energies_found++; // we found a new energy!
      if (max_energy_seen < 0 || E > max_energy_seen) max_energy_seen = E;
      if (min_energy_seen < 0 || E < min_energy_seen) min_energy_seen = E;
      if (use_sad) {
        printf("  (moves %ld, energies_found %ld, erange: %d -> %d effective t0 = %g)\n",
               moves, energies_found, min_energy_seen, max_energy_seen,
               sa_prefactor*energies_found
                 *(max_energy_seen-min_energy_seen)/(minT*use_sad));
      }
    }
    if (use_sad && energies_found > 1) {
      wl_factor = sa_prefactor*energies_found*(max_energy_seen-min_energy_seen)
        /(minT*moves*use_sad);
    } else {
      wl_factor = sa_prefactor*sa_t0/max(sa_t0, moves);
    }
  }
  energy_histogram[E]++;
  //if(pessimistic_observation[min_important_energy]) walkers_up[E]++;
  if (use_sad && energies_found > 1 && wl_factor > 0) {
    if (E < too_high_energy) {
      // We are at higher energy than the maximum entropy state, so we
      // need to tweak our weights by even more, since we don't spend
      // much time here.

      ln_energy_weights[E] =
        min(ln_energy_weights[E] - wl_factor,
            ln_energy_weights[too_high_energy] -
            log(wl_factor
                + exp(ln_energy_weights[too_high_energy]-ln_energy_weights[E])));
      if (!(isnormal(ln_energy_weights[E]) || ln_energy_weights[E] == 0)) {
        printf("lnw[%d] = %g\n", E, ln_energy_weights[E]);
      }
      assert(isnormal(ln_energy_weights[E]) || ln_energy_weights[E] == 0);
    } else if (E > too_low_energy) {

      ln_energy_weights[E] =
        min(ln_energy_weights[E] - wl_factor,
            ln_energy_weights[too_low_energy] + (E-too_low_energy)/minT
            - log(wl_factor + exp(ln_energy_weights[too_low_energy] - ln_energy_weights[E] + (E-too_low_energy)/minT)));

      if (ln_energy_weights[E] < ln_energy_weights[min_important_energy] + (E-min_important_energy)/minT) {
        ln_energy_weights[E] =
          ln_energy_weights[min_important_energy] + (E-min_important_energy)/minT;
        too_low_energy = E;
        // printf("We were almost very cray at energy %d\n", energy);
      }

      if (!(isnormal(ln_energy_weights[E]) || ln_energy_weights[E] == 0)) {
        printf("lnw[%d] = %g\n", E, ln_energy_weights[E]);
      }
      assert(isnormal(ln_energy_weights[E]) || ln_energy_weights[E] == 0);
    } else {
      // We are in the "interesting" region, so use an ordinary SA update.
      ln_energy_weights[E] -= wl_factor;
    }
  } else {
    // Note: if not using WL or SA method, wl_factor = 0 and the
    // following has no effect.
    ln_energy_weights[E] -= wl_factor;
  }
  if (use_sad) {
    bool print_edges = false;
    if (ln_energy_weights[E] < ln_energy_weights[max_entropy_state]
        || energy_histogram[max_entropy_state] == 0) {
      // We are now at the max_entropy_state!
      max_entropy_state = E;
      if (max_entropy_state < too_high_energy) too_high_energy = max_entropy_state;
      // print_edges = true;
    }
    if (ln_energy_weights[E] < ln_energy_weights[min_important_energy] + (E-min_important_energy)/minT
        || energy_histogram[min_important_energy] == 0) {
      // We are above the min_important_energy tangent line, which
      // means we are the new min_important_energy.
      //if (energy != min_important_energy) print_edges = true;
      min_important_energy = E;
      if (min_important_energy > too_low_energy) too_low_energy = min_important_energy;
    }

    if (print_edges) printf(" %5d ...%5d -->%5d ...%5d\n",
                            too_low_energy, min_important_energy, max_entropy_state, too_high_energy);
  }
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

void ising_simulation::compute_ln_dos(dos_types dos_type) {
  if (dos_type == histogram_dos) {
    for (int i = max_entropy_state; i < energy_levels; i++) {
      if (energy_histogram[i] != 0) {
        ln_dos[i] = log(energy_histogram[i]) + ln_energy_weights[i];
      } else {
        ln_dos[i] = -DBL_MAX; // located in <float.h>.
      }
    }
  } else if (dos_type == weights_dos) {
    // weights_dos is useful for WL, SAD, or SAMC algorithms, where
    // the density of states is determined directly from the weights.
    int minE = 0;
    double betamax = 1.0/minT;
    if (use_sad) {
      // This is hokey, but we just want to ensure that we have a proper
      // max_entropy_state before computing the ln_dos.
      max_entropy_state = 0;
      for (int i=0; i<energy_levels; i++) {
        if (ln_energy_weights[i] < ln_energy_weights[max_entropy_state]) {
          max_entropy_state = i;
        }
      }
    }
    for (int i=0; i<energy_levels; i++) {
      ln_dos[i] = ln_energy_weights[max_entropy_state] - ln_energy_weights[i];
      if (!minE && ln_dos[i-1] - ln_dos[i] > betamax) minE = i;
    }
    if (use_sad) {
      // Above the max_entropy_state our weights are effectively constant,
      // so the density of states is proportional to our histogram.
      //~ for (int i=0; i<max_entropy_state; i++) {
        //~ if (energy_histogram[i]) {
          //~ ln_dos[i] = log(energy_histogram[i]/double(energy_histogram[max_entropy_state]));
        //~ }
      //~ }
      // Below the minimum important energy, we also need to use the histogram,
      // only now adjusted by a Boltzmann factor.  We compute the min
      // important energy above for extreme clarity.
      //~ for (int i=minE+1; i<energy_levels; i++) {
        //~ if (energy_histogram[i]) {
          //~ ln_dos[i] = ln_dos[minE] + log(energy_histogram[i]/double(energy_histogram[minE]))
            //~ - (i-minE)*betamax;  // the last bit gives Boltzmann factor
        //~ }
      //~ }
      // Now let us set the ln_dos for any sites we have never visited
      // to be equal to the minimum value of ln_dos for sites we
      // *have* visited.  These are unknown densities of states, and
      // there is no particular reason to set them to be crazy low.
      double lowest_ln_dos = 0;
      for (int i=0; i<energy_levels; i++) {
        if (energy_histogram[i]) lowest_ln_dos = min(lowest_ln_dos, ln_dos[i]);
      }
      for (int i=0; i<energy_levels; i++) {
        if (!energy_histogram[i]) ln_dos[i] = lowest_ln_dos;
      }
    }
  }
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

  double fix_kT = 0;
  int samc = false;
  int sad = false;

  int NN = 10;
  int resume = false;
  long total_moves = 10000;

  char *data_dir = new char[1024];
  sprintf(data_dir,"none");
  const char *default_data_dir = "papers/histogram/data/ising";
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

    {"sad", '\0', POPT_ARG_NONE, &sad, 0,
     "Use stochastic approximation monte carlo dynamical version", "BOOLEAN"},

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

  // Check that only one histogram method is used
  if (sad + samc + (fix_kT != 0) + resume != 1) {
    printf("Exactly one histogram method must be selected! (%d %d %g)\n",
           sad, samc, fix_kT);
    return 254;
  }
  // Set default data directory
  if (strcmp(data_dir,"none") == 0) {
    sprintf(data_dir,"%s",default_data_dir);
    printf("\nUsing default data directory: [deft]/%s\n",data_dir);
  }

  mkdir(data_dir, 0777); // create save directory

  char *ising_fname = new char[1024];
  sprintf(ising_fname, "%s/%s.dat", data_dir, filename);

  ising_simulation ising(NN);
  ising.canonical_temperature(6);

  took("Starting program");

  // ----------------------------------------------------------------------------
  // Resume functionality
  // ----------------------------------------------------------------------------

  if (resume) {
     // We are continuing a previous simulation.
    if (sad) {// This will be a Monte Carlo method eventually.
      // Open the file!
      FILE *rfile = fopen((const char *)ising_fname, "r");
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
        int resumeNN = 0;
        if (fscanf(rfile, " NN = %d\n", &resumeNN) != 1) {
          printf("error reading NN!\n");
          exit(1);
        }
        if (resumeNN != NN) {
          printf("Error:  must specify same N value when resuming! %d != %d\n",
                 NN, resumeNN);
          exit(1);
        }
        printf("NN is now %d\n", NN);
        if (fscanf(rfile, "total-moves = %ld\n", &total_moves) != 1) {
          printf("error reading total-moves!\n");
          exit(1);
        }
        printf("total-moves is now %ld\n", total_moves);
        if (fscanf(rfile, "minT = %lg\n", &minT) != 1) {
          printf("error reading minT!\n");
          exit(1);
        }
        if (fscanf(rfile, "moves = %ld\n", &ising.moves) != 1) {
          printf("error reading ising.moves!\n");
          exit(1);
        }
        if (fscanf(rfile, "E = %d\n", &ising.E) != 1) {
          printf("error reading E!\n");
          exit(1);
        }
        if (fscanf(rfile, "J = %d\n", &J) != 1) {
          printf("error reading J!\n");
          exit(1);
        }

        printf("lndos is now {\n");
        fscanf(rfile, " lndos = np.array([");
        for (int i = 0; i < ising.energy_levels; i++) {
          if (fscanf(rfile,"\t%lg,\n", &ising.ln_dos[i]) != 1) {
            printf("error reading lndos at energy #%d!\n", i);
            exit(1);
          }
        }
        fscanf(rfile, " ])\n");
        printf("}\n");

        fscanf(rfile, " lnw = np.array([");
        for (int i = 0; i < ising.energy_levels; i++) {
          if (fscanf(rfile,"\t%lg,\n", &ising.ln_energy_weights[i]) != 1) {
            printf("error reading lnw at energy #%d!\n", i);
            exit(1);
          }
        }
        fscanf(rfile, " ])\n");

        fscanf(rfile, " histogram = np.array([");
        for (int i = 0; i < ising.energy_levels; i++) {
          if (fscanf(rfile,"\t%ld,\n", &ising.energy_histogram[i]) != 1) {
            printf("error reading histogram at energy #%d!\n", i);
            exit(1);
          }
        }
        fscanf(rfile, " ])\n");

        fscanf(rfile, " S = np.array([\n");
        for (int i=0; i<NN; i++) {
          fscanf(rfile, "\t[");
          for (int j=0; j<NN; j++) {
            fscanf(rfile,"%2d,", &ising.S[i+NN*j]);
          }
          fscanf(rfile, "],\n");
        }
        fscanf(rfile, "])\n");

        fclose(rfile);
      }
      } else {
        printf("I do not know how to resume yet!\n");
        exit(1);
      }
      took("Reading resume file");
  }

  printf("version: %s\n",version_identifier());

  ising.calculate_energy();
  for (long i = 0; i < total_moves; i++) {
    ising.flip_a_spin();
  }

  printf("I think the energy is %d\n", ising.E);
  ising.calculate_energy();
  printf("the energy is %d\n", ising.E);

  took("Running");

  ising.compute_ln_dos(histogram_dos);

  for(int i = 0; i < ising.energy_levels; i++){
    if (ising.energy_histogram[i] > 0) {
      ising.energies_found++;
      printf("energy histogram at %ld is %ld with lndos %g\n",
              ising.energies_found,ising.energy_histogram[i],ising.ln_dos[i]);
    }
  }

  // ----------------------------------------------------------------------------
  // Generate save file info
  // ----------------------------------------------------------------------------

  // Save resume file
      {
        FILE *ising_out = fopen((const char *)ising_fname, "w");

        fprintf(ising_out, "import numpy as np\n\n");
        random::dump_resume_info(ising_out);
        fprintf(ising_out,"version = %s\n\n",version_identifier());

        if (fix_kT) {
          fprintf(ising_out,"method = 'canonical'\nkT = %g\n", fix_kT);
        } else if (samc) {
            fprintf(ising_out,"method = 'samc'\n");
            // fprintf(ising_out,"samc_t0 = %g\n",ising.sa_t0);
        } else if (sad) {
            fprintf(ising_out,"method = 'sad'\n");
            // fprintf(ising_out,"too_high_energy = %d\n", ising.too_high_energy);
        }

        fprintf(ising_out,"NN = %i\n",NN);
        fprintf(ising_out,"total-moves = %li\n",total_moves);
        fprintf(ising_out,"minT = %g\n",minT);
        fprintf(ising_out,"moves = %li\n",ising.moves);
        fprintf(ising_out,"E = %i\n",ising.E);
        fprintf(ising_out,"J = %i\n", J);
        // inserting arrays into text file.
        fprintf(ising_out, "lndos = np.array([\n");
        for (int i = 0; i < ising.energy_levels; i++) {
          fprintf(ising_out,"\t%.17g,\n", ising.ln_dos[i]);
        }
        fprintf(ising_out, "])\n");
        fprintf(ising_out, "lnw = np.array([\n");
        for (int i = 0; i < ising.energy_levels; i++) {
          fprintf(ising_out,"\t%.17g,\n", ising.ln_energy_weights[i]);
        }
        fprintf(ising_out, "])\n");
        fprintf(ising_out, "histogram = np.array([\n");
        for (int i = 0; i < ising.energy_levels; i++) {
          fprintf(ising_out,"\t%ld,\n", ising.energy_histogram[i]);
        }
        fprintf(ising_out, "])\n");
        fprintf(ising_out, "S = np.array([\n");
        for (int i=0; i<NN; i++) {
          fprintf(ising_out, "\t[");
          for (int j=0; j<NN; j++) {
            fprintf(ising_out,"%2d,", ising.S[i+NN*j]);
          }
          fprintf(ising_out, "],\n");
        }
        fprintf(ising_out, "])\n");

        fclose(ising_out);
      }

  // -------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // -------------------------------------------------------------------
}

