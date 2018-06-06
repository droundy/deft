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

// TODO items:

// 1. Write a script that tests whether resume works.  (And fix
//    resume.)  i.e. the script will run the code once with fixed #
//    moves, then run again resuming.  Then it'll run with twice the
//    moves and confirm the results agree.  Eventually we'll need the
//    test to work with all methods.

// 2. Make sure the min_important_energy is always up-to-date (but
//    keep it cheap).

// 3. Make sure the max_entropy_energy is always up-to-date (but keep
//    it cheap).

// 4. Fix FIXMEs.

// ---------------------------------------------------------------------
// Define "Constants" -- set from arguments then unchanged
// ---------------------------------------------------------------------

const int Q = 2;   // Define number of spin states.

// ---------------------------------------------------------------------
// ising.h files -- this could likely go into an ising.h file
// ---------------------------------------------------------------------

struct energy {
  int value;
  explicit energy(int val) : value(val) {}
  bool operator<(energy other) const { return value < other.value; }
  bool operator<=(energy other) const { return value <= other.value; }
  bool operator>(energy other) const { return value > other.value; }
  bool operator>=(energy other) const { return value >= other.value; }
  void operator+=(energy other) { value += other.value; }
  energy operator+(energy other) const { return energy(value + other.value); }
  void operator-=(energy other) { value -= other.value; }
  energy operator-(energy other) const { return energy(value - other.value); }
};

enum dos_types { histogram_dos, transition_dos, weights_dos };

struct simulation_parameters {
  int N;            // N*N is the number of sites

  int J;       // Define our spin coupling.
  double minT; // Define minimum Temperature.

  /* The following determine how we perform an ising flip */
  double sa_t0; // if it is non-zero, update gamma on each move SA-style.
  double sa_prefactor; // prefactor in computing gamma when running SA.

  int use_sad;  // use sad method.
  int use_wl; // if true, we are using WL

  double T;      // canonical temperature.
  simulation_parameters() {
    N=0;
    J=1;
    T = 0;
    minT=0.2;
    sa_t0 = 0;
    sa_prefactor = 1;
    use_sad = false;
    use_wl = false;
  }
};


struct ising_simulation {
  simulation_parameters param; // the parameters that specify how to simulate

  long moves;   // the current number of moves
  int *S;
  energy E;      // system energy.

  // the last time we printed status text (i.e. from initialization)
  double estimated_time_per_iteration; // in units of seconds per iteration

  /* The following determine how we perform an ising flip */
  double gamma; // if it is non-zero, update weights on each move WL-style.

  long energy_levels;

  long energies_found;
  long *energy_histogram;

  double *ln_energy_weights;
  double *ln_dos;
  energy max_entropy_energy, min_important_energy;
  energy min_energy; // The lowest energy we have ever found for this
                  // system.  This must be kept up-to-date because
                  // it is used by SAD to compute the delta E.
  energy max_energy; // The highest energy we have ever found for this
                  // system.  This must be kept up-to-date because
                  // it is used by SAD to compute the delta E.
  energy too_hi_energy; // too_hi_energy is used in SAD, and is the
                       // highest ever value that max_entropy_energy has
                       // taken.  It is used to define the range over
                       // which our histogram is made flat.
  energy too_lo_energy; // too_lo_energy is used in SAD, and is the
                      // lowest ever value that min_important_energy has
                      // taken.  It is used to define the range over
                      // which our histogram is made flat.

  explicit ising_simulation(simulation_parameters par); // generate the spin lattice
  ~ising_simulation();

  int random_flip(int oldspin) const; // pick a new random (but changed) spin

  void flip_a_spin();
  void end_flip_updates();
  double calculate_energy();

  long index_from_energy(energy E) const {
    return param.J*param.N*param.N/2 + E.value/4;
  }
  energy energy_from_index(int i) const {
    return energy(4*(i - param.J*param.N*param.N/2));
  }

  void compute_ln_dos(dos_types dos_type);
  void initialize_samc(int am_sad);
  bool reached_iteration_cap();
};

// ising_simulation methods

int ising_simulation::random_flip(int oldspin) const {
  if (Q==2) return -oldspin;

  return random::ran64() % Q - Q/2; // gives value -Q/2 to Q/2-1?
}

ising_simulation::ising_simulation(simulation_parameters par)
  : param(par), E(0), max_entropy_energy(0), min_important_energy(0),
    min_energy(0), max_energy(0), too_hi_energy(0), too_lo_energy(0) {

  // seconds per iteration (will be adjusted from actual timing)
  estimated_time_per_iteration = 0.1;
  moves = 0;
  const int J = param.J;
  const int N = param.N;
  S = new int[N*N];
  // energy histogram
  energy_levels = J*N*N;
  ln_energy_weights = new double[energy_levels]();
  ln_dos = new double[energy_levels]();
  energy_histogram = new long[energy_levels]();
  //printf("energy_levels %ld\n", energy_levels);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      S[i+N*j] = 1; // initialize to all 1s
    }
  }
  calculate_energy();
  min_energy = E;
  max_energy = E;
  max_entropy_energy = E;
  too_hi_energy = E;
  too_lo_energy = E;
  energies_found = 1; // we found just one energy

  if (param.T > 0) {
    // We are doing a canonical simulation, so set up the weights
    // appropriately!
    for (int i=0; i<energy_levels; i++) {
      ln_energy_weights[i] = 4*i/param.T;
    }
  }

}

ising_simulation::~ising_simulation() {
  delete[] energy_histogram;
  delete[] ln_energy_weights;
  delete[] ln_dos;
  delete[] S;
}

void ising_simulation::flip_a_spin() {
  const int J = param.J;
  const int N = param.N;
  moves += 1;
  int k = random::ran64() % (N*N);
  int i = k % N;
  int j = k / N;
  int old = S[j + i*N];
  S[j + i*N] = random_flip(S[j + i*N]);

  int neighbor_spins = S[j + ((i+1) % N)*N] + S[j + ((i-1+N) % N)*N] +
                       S[((j+1) % N) + i*N] + S[((j-1+N) % N) + i*N];

  const energy deltaE = energy(J*(old - S[j + i*N])*neighbor_spins);

  double lnprob = 0;
  if (param.use_sad) {
    const energy e1 = E;
    const energy e2 = E + deltaE;
    double lnw1 = ln_energy_weights[index_from_energy(e1)];
    if (e1 > too_hi_energy) {
      lnw1 = ln_energy_weights[index_from_energy(too_hi_energy)];
    } else if (e1 < too_lo_energy) {
      lnw1 = ln_energy_weights[index_from_energy(too_lo_energy)]
             - (too_lo_energy - e1).value/param.minT;
    }
    double lnw2 = ln_energy_weights[index_from_energy(e2)];
    if (e2 > too_hi_energy) {
      lnw2 = ln_energy_weights[index_from_energy(too_hi_energy)];
    } else if (e2 < too_lo_energy) {
      lnw2 = ln_energy_weights[index_from_energy(too_lo_energy)]
             - (too_lo_energy - e2).value/param.minT;
    }
    lnprob = lnw1 - lnw2;
  } else {
    lnprob = ln_energy_weights[index_from_energy(E)]
              - ln_energy_weights[index_from_energy(E+deltaE)];
  }
  if (lnprob < 0 && random::ran() > exp(lnprob)) {
    S[j + i*N] = old;
    //printf("not flipping from E=%g\n", E);
  } else {
    E += deltaE;
  }

  energy_histogram[index_from_energy(E)] += 1;
  end_flip_updates();

}

void ising_simulation::end_flip_updates(){
  // update iteration counter, energy histogram, etc
  if (param.sa_t0 || param.use_sad) {
    if (energy_histogram[index_from_energy(E)] == 0) {
      energies_found++; // we found a new energy!
      if (E > max_energy) max_energy = E;
      if (E < min_energy) min_energy = E;
      if (param.use_sad) {
        printf("  (moves %ld, energies_found %ld, erange: %d -> %d)\n",
               moves, energies_found, min_energy.value, max_energy.value);
      }
    }
    if (param.use_sad && energies_found > 1) {
      gamma = param.sa_prefactor*energies_found*(max_energy.value-min_energy.value)
        /(3*param.minT*moves);
    } else {
      gamma = param.sa_prefactor*param.sa_t0/max(param.sa_t0, moves);
    }
  }
  energy_histogram[index_from_energy(E)]++;
  if (param.use_sad && energies_found > 1) {
    if (E > too_hi_energy) {
      // We are at higher energy than the maximum entropy state, so we
      // need to tweak our weights by even more, since we don't spend
      // much time here.
      ln_energy_weights[index_from_energy(E)] =
        max(ln_energy_weights[index_from_energy(E)] + gamma,
            ln_energy_weights[index_from_energy(too_hi_energy)] +
            log(gamma
                + exp(ln_energy_weights[index_from_energy(E)]
                      - ln_energy_weights[index_from_energy(too_hi_energy)])));
      if (!(isnormal(ln_energy_weights[index_from_energy(E)])
            || ln_energy_weights[index_from_energy(E)] == 0)) {
        printf("lnw[%d] = %g\n", E.value, ln_energy_weights[index_from_energy(E)]);
      }
      assert(isnormal(ln_energy_weights[index_from_energy(E)])
              || ln_energy_weights[index_from_energy(E)] == 0);
    } else if (E < too_lo_energy) {

      ln_energy_weights[index_from_energy(E)] =
        max(ln_energy_weights[index_from_energy(E)] + gamma,
            ln_energy_weights[index_from_energy(too_lo_energy)]
            + (too_lo_energy - E).value/param.minT
            + log(gamma + exp(ln_energy_weights[index_from_energy(E)]
                              - ln_energy_weights[index_from_energy(too_lo_energy)]
                              + (too_lo_energy - E).value/param.minT)));

      if (ln_energy_weights[index_from_energy(E)] >
            ln_energy_weights[index_from_energy(min_important_energy)]
            - (min_important_energy - E).value/param.minT) {
        // FIXME Think about whether this is needed.
        ln_energy_weights[index_from_energy(E)] =
            ln_energy_weights[index_from_energy(min_important_energy)]
            - (min_important_energy - E).value/param.minT;
        too_lo_energy = E;
        // printf("We were almost very cray at energy %d\n", energy);
      }

      if (!(isnormal(ln_energy_weights[index_from_energy(E)])
            || ln_energy_weights[index_from_energy(E)] == 0)) {
        printf("lnw[%d] = %g\n", E.value, ln_energy_weights[index_from_energy(E)]);
      }
      assert(isnormal(ln_energy_weights[index_from_energy(E)])
              || ln_energy_weights[index_from_energy(E)] == 0);
    } else {
      // We are in the "interesting" region, so use an ordinary SA update.
      ln_energy_weights[index_from_energy(E)] += gamma;
    }
  } else {
    // Note: if not using WL or SA method, gamma = 0 and the
    // following has no effect.
    ln_energy_weights[index_from_energy(E)] += gamma;
  }
  if (param.use_sad) {
    bool print_edges = false;
    if (ln_energy_weights[index_from_energy(E)] >
        ln_energy_weights[index_from_energy(max_entropy_energy)]
        || energy_histogram[index_from_energy(max_entropy_energy)] == 0) {
      // We are now at the max_entropy_energy!
      max_entropy_energy = E;
      if (max_entropy_energy > too_hi_energy) too_hi_energy = max_entropy_energy;
      // print_edges = true;
    }
    if (ln_energy_weights[index_from_energy(E)] >
        ln_energy_weights[index_from_energy(min_important_energy)]
        - (min_important_energy - E).value/param.minT
        || energy_histogram[index_from_energy(min_important_energy)] == 0) {
      // We are above the min_important_energy tangent line, which
      // means we are the new min_important_energy.
      //if (energy != min_important_energy) print_edges = true;
      min_important_energy = E;
      if (min_important_energy < too_lo_energy) too_lo_energy = min_important_energy;
    }

    if (print_edges) printf(" %5d ...%5d -->%5d ...%5d\n",
                            too_lo_energy.value, min_important_energy.value,
                            max_entropy_energy.value, too_hi_energy.value);
  }
}

double ising_simulation::calculate_energy() {
  const int N = param.N;
  const int J = param.J;
  E.value = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      // we want to enforce periodic boundary conditions.
      int neighbor_spins = S[(i+1) % N + j*N] + S[i + ((j+1) % N)*N];

      E.value += -J*neighbor_spins*S[i+j*N];
    }
  }
  return E.value;
}

void ising_simulation::compute_ln_dos(dos_types dos_type) {
  if (dos_type == histogram_dos) {
    for (int i = index_from_energy(max_entropy_energy); i < energy_levels; i++) {
      if (energy_histogram[i] != 0) {
        ln_dos[i] = log(energy_histogram[i]) + ln_energy_weights[i];
      } else {
        ln_dos[i] = -DBL_MAX; // located in <float.h>.
      }
    }
  } else if (dos_type == weights_dos) {
    // weights_dos is useful for WL, SAD, or SAMC algorithms, where
    // the density of states is determined directly from the weights.
    double max_entropy = ln_energy_weights[index_from_energy(max_entropy_energy)];
    for (int i=0; i<energy_levels; i++) {
      //if (energy_histogram[i] != 0) {FIXME: WORKING HERE NOW!!!
      ln_dos[i] = ln_energy_weights[i] - max_entropy;
    }
    if (param.use_sad) {
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

  simulation_parameters param;

  int resume = false;
  long total_moves = 10000;

  char *data_dir = new char[1024];
  sprintf(data_dir,"none");
  const char *default_data_dir = "papers/histogram/data/ising";
  char *filename = new char[1024];
  sprintf(filename, "none");
  char *filename_suffix = new char[1024];
  sprintf(filename_suffix, "none");

  poptContext optCon;

  // -------------------------------------------------------------------
  // Parse input options
  // -------------------------------------------------------------------

  poptOption optionsTable[] = {
    //{"seed", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0,
    // "Seed for the random number generator", "INT"},

    /*** ISING MODEL PARAMETERS ***/

    {"N", '\0', POPT_ARG_INT, &param.N, 0, "N*N is the number of spin sites", "INT"},

    /*** SIMULATION ITERATIONS ***/

    {"total-moves", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &total_moves,
     0, "Number of moves for which to run the simulation", "INT"},

    /*** PARAMETERS DETERMINING OUTPUT FILE DIRECTORY AND NAMES ***/

    {"dir", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &data_dir, 0,
     "Directory in which to save data", "data_dir"},
    {"filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     "Base of output file names", "STRING"},


    /*** HISTOGRAM METHOD OPTIONS ***/
    {"resume", '\0', POPT_ARG_NONE, &resume, 0,
     "Resume previous simulation", "BOOLEAN"},

    {"sad", '\0', POPT_ARG_NONE, &param.use_sad, 0,
     "Use stochastic approximation monte carlo dynamical version", "BOOLEAN"},
    {"wl", '\0', POPT_ARG_NONE, &param.use_wl, 0,
     "Use Wang-Landau method", "BOOLEAN"},
    {"minT", '\0', POPT_ARG_DOUBLE, &param.minT, 0,
     "The minimum temperature we care about", "DOUBLE"},
    {"sa-t0", '\0', POPT_ARG_DOUBLE,
     &param.sa_t0, 0, "use SAMC with specified t0 value", "DOUBLE"},

    /*** HISTOGRAM METHOD PARAMETERS ***/

    //{"gamma", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gamma,
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
  if (param.use_sad + (param.sa_t0 != 0) + (param.T != 0) != 1) {
    printf("Exactly one histogram method must be selected! (%d %g %g)\n",
           param.use_sad, param.sa_t0, param.T);
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

  ising_simulation ising(param);
  took("Starting program");

  // ----------------------------------------------------------------------------
  // Resume functionality
  // ----------------------------------------------------------------------------

  if (resume) {
    // We are continuing a previous simulation.
    FILE *rfile = fopen((const char *)ising_fname, "r");
    if (param.use_sad) {
      // Open the file!
      //FILE *rfile = fopen((const char *)ising_fname, "r");
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
        {
          int resumeN = 0;
          if (fscanf(rfile, " N = %d\n", &resumeN) != 1) {
            printf("error reading N!\n");
            exit(1);
          }
          if (resumeN != param.N) {
            printf("Error:  must specify same N value when resuming! %d != %d\n",
                   param.N, resumeN);
            exit(1);
          }
        }
        if (fscanf(rfile, "total-moves = %ld\n", &total_moves) != 1) {
          printf("error reading total-moves!\n");
          exit(1);
        }
        if (fscanf(rfile, "moves = %ld\n", &ising.moves) != 1) {
          printf("error reading ising.moves!\n");
          exit(1);
        }
        if (fscanf(rfile, "E = %d\n", &ising.E.value) != 1) {
          printf("error reading E!\n");
          exit(1);
        }
        {
          int resumeJ;
          if (fscanf(rfile, "J = %d\n", &resumeJ) != 1) {
            printf("error reading J!\n");
            exit(1);
          }
          if (resumeJ != param.J) {
            printf("Error:  must specify same J value when resuming! %d != %d\n",
                   param.J, resumeJ);
            exit(1);
          }
        }
        if (fscanf(rfile, "E_found = %li\n", &ising.energies_found) != 1) {
          printf("error reading E_found!\n");
          exit(1);
        }
        if (param.use_sad) {
          fscanf(rfile," too_hi_energy = %i\n", &ising.too_hi_energy.value); // NAME?
          fscanf(rfile," too_lo_energy = %i\n", &ising.too_lo_energy.value); // NAME?
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
        for (int i=0; i<param.N; i++) {
          fscanf(rfile, "\t[");
          for (int j=0; j<param.N; j++) {
            fscanf(rfile,"%2d,", &ising.S[i+param.N*j]);
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

  // MAIN CODE EXECUTION HERE!
  ising.calculate_energy();
  for (long i = 0; i < total_moves; i++) {
    ising.flip_a_spin();
  }

  printf("I think the energy is %d\n", ising.E.value);
  ising.calculate_energy();
  printf("the energy is %d\n", ising.E.value);

  took("Running");

  if (param.use_sad) {
    ising.compute_ln_dos(weights_dos);
  } else {
    ising.compute_ln_dos(histogram_dos);
  }

  for(int i = 0; i < ising.energy_levels; i++){
    if (ising.energy_histogram[i] > 0) {
      ising.energies_found++;
      if (ising.energy_from_index(i) > ising.max_energy) {
       ising.max_energy = ising.energy_from_index(i);
      }
      if (ising.energy_from_index(i) < ising.min_energy) {
       ising.min_energy = ising.energy_from_index(i);
      }
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

    if (param.T) {
      fprintf(ising_out,"method = 'canonical'\nkT = %g\n", param.T);
    } else if (param.sa_t0) {
      fprintf(ising_out,"method = 'samc'\n");
      // fprintf(ising_out,"samc_t0 = %g\n",ising.sa_t0);
    } else if (param.use_sad) {
      fprintf(ising_out,"method = 'sad'\n");
    }

    fprintf(ising_out,"N = %i\n",param.N);
    fprintf(ising_out,"minT = %g\n",param.minT);
    fprintf(ising_out,"moves = %li\n",ising.moves);
    fprintf(ising_out,"E = %i\n",ising.E.value);
    fprintf(ising_out,"J = %i\n", param.J);
    fprintf(ising_out,"E_found = %li\n", ising.energies_found); // NAME?
    if (param.use_sad) {
      fprintf(ising_out,"too_hi_energy = %i\n", ising.too_hi_energy.value); // NAME?
      fprintf(ising_out,"too_lo_energy = %i\n", ising.too_lo_energy.value); // NAME?
    }
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
    for (int i=0; i<param.N; i++) {
      fprintf(ising_out, "\t[");
      for (int j=0; j<param.N; j++) {
        fprintf(ising_out,"%2d,", ising.S[i+param.N*j]);
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

