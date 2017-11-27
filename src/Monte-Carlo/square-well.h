#include "../vector3d.h"
#pragma once

//#define NDEBUG // disable assertions for extra speed! (and to test if it matters)

struct ball {
  vector3d pos;
  double R;
  int *neighbors;
  int num_neighbors;
  vector3d neighbor_center;

  ball();
  ball(const ball &p);

  ball operator=(const ball &p);
};

// Struct to store all information about an appempt to move one ball
struct move_info {
  long total_old;
  long total;
  long working_old;
  long working;
  int updates;
  int informs;
  move_info();
};

enum end_conditions { none, optimistic_min_samples, pessimistic_min_samples,
                      optimistic_sample_error, pessimistic_sample_error, flat_histogram,
                      init_iter_limit };

enum dos_types { histogram_dos, transition_dos, weights_dos };

// This should store all information needed to run a simulation.  Thus
// we can just pass this struct around to functions that run the
// simulation.  We do not maintain here any of the histograms except
// the energy_histogram, since only that histogram is needed to
// actually run the simulation.  Collecting "real" data is another
// matter.
struct sw_simulation {
  long iteration; // the current iteration number

  /* The following describe the current state of the system. */

  ball *balls;
  int energy;

  /* The following are constant parameters that describe the physical
     system, but do not change as we simulate. */
  double well_width; // width of well in units of sphere radius
  double filling_fraction; // proportion of space filled with spheres
  int N; // number of balls
  double len[3]; // the size of the cell
  int walls; // should this be an enum?
  int sticky_wall; // do we have an attractive region near one wall?
  double interaction_distance; // the distance over which balls interact

  /* The following are constant parameters that describe the details
     of how we run the computation, and do not change as we
     simulate. */

  // max_neighbors is the upper limit to the maximum number of
  // neighbors a ball could have.  It is needed as an imput to
  // move_a_ball.
  int max_neighbors;
  double neighbor_R; // radius of our neighbor sphere
  double translation_scale; // scale for how far to move balls
  int energy_levels; // total number of energy levels
  int energies_found; // number of energy levels that have been visited

  /* The following accumulate results of the simulation. Although
     ln_energy_weights is a constant except during initialization. */

  int max_entropy_state, min_energy_state, min_important_energy;
  move_info moves;
  long *energy_histogram;
  double *ln_energy_weights;
  dos_types sim_dos_type;

  /* The following keep track of how many times we have walked
     between the a given energy and the state of max entropy */

  long *optimistic_samples; // how many times have we gone down to this energy?
  long *pessimistic_samples; // how many samples from the point of max entropy have we had?
  bool *pessimistic_observation; // of a given energy

  /* The following control end conditions for histogram methods */

  end_conditions end_condition;
  double max_time, start_time; // maximum time to run simulation, and when we started
  double min_T; // minimum temperature we care about
  int min_samples; // force some number of minimum energy samples
  double sample_error; // the maximum fractional sample error to achieve in initialization
  double flatness; // maximum allowable proportional deviation from mean histogram value
  int init_iters; // number of iterations for which to initialize

  /* The following determine how we move balls */
  bool use_tmmc; // true means rejection/acceptance comes from transition matrix.
  double wl_factor; // if it is non-zero, update weights on each move WL-style.
  double sa_t0; // if it is non-zero, update wl_factor on each move SA-style.
  double sa_prefactor; // prefactor in computing wl_factor when running SA.
  bool use_satmmc; // if true, reset sa_t0 whenever we encounter a new energy (also do TMMC).
  bool use_sad;  // if true, dynamically update sa_t0 when we encounter new energies.
  double sa_weight; // the most recent sa_weight used in satmmc.

  /* The following define file names for periodic output files that
     are dumped every so often.  It should contain a single %d style format. */
  const char *transitions_movie_filename_format;
  mutable int transitions_movie_count;
  const char *dos_movie_filename_format;
  mutable int dos_movie_count;
  const char *lnw_movie_filename_format;
  mutable int lnw_movie_count;
  /* Finally, the following defines file names for output files. */
  const char *transitions_filename;

  /* The following tracks how many transitions we have attempted from
     a given energy level to nearby energy levels.  The advantage of
     this metric is that we can (hopefully) get away without resetting
     it when we change the weights.  Thus we could accumulate better
     statistics on entropy differences, under the assumption that we
     sample all states of a given energy equally. */
  int biggest_energy_transition;
  long *collection_matrix;
  long &transitions(int energy, int energy_change) {
    assert(energy_change >= -biggest_energy_transition);
    assert(energy_change <= biggest_energy_transition);
    assert(energy >= 0);
    assert(energy < energy_levels);
    return collection_matrix[energy*(2*biggest_energy_transition+1)
                             + energy_change+biggest_energy_transition];
  };
  long transitions(int energy, int energy_change) const {
    assert(energy_change >= -biggest_energy_transition);
    assert(energy_change <= biggest_energy_transition);
    assert(energy >= 0);
    assert(energy < energy_levels);
    return collection_matrix[energy*(2*biggest_energy_transition+1)
                             + energy_change+biggest_energy_transition];
  };
  /* "transition_matrix" is a read-only sloppy and normalized version
     of the matrix also called "transitions" above, which is a little
     easier for me to wrap my brains around.  DJR */
  double transition_matrix(int to, int from) const {
    if (abs(to - from) > biggest_energy_transition ||
        to < 0 || from < 0 || to >= energy_levels || from >= energy_levels) {
      return 0;
    }
    long norm = 0;
    for (int de=-biggest_energy_transition; de<=biggest_energy_transition; de++) {
      norm += transitions(from, de);
    }
    if (norm == 0) return 0;
    return transitions(from, to - from)/double(norm);
  };

  int max_interactions() const { // return the maximum observed number of interactions
    for (int i = energy_levels-1; i >= 0; i--) if (energy_histogram[i]) return i;
    return 0;
  };

  /* Up-moving walkers for optimized ensemble method */
  long *walkers_up;

  void reset_histograms();
  void move_a_ball(); // attempt to move one ball
  void end_move_updates(); // updates to run at the end of every move
  void energy_change_updates(int energy_change); // updates to run if we've changed energy

  // the last time we printed status text (i.e. from initialization)
  double estimated_time_per_iteration; // in units of seconds per iteration

  // iterate long enough to find the max entropy state and initialize
  // the translation distance. return most probable energy
  int initialize_max_entropy(double acceptance_goal = 0.4);

  // initialize the translation distance. return most probable energy
  void initialize_translation_distance(double acceptance_goal = 0.4);

  // iterate enough times for the energy to change n times.  Return
  // the number of "up" moves.
  int simulate_energy_changes(int num_moves);

  // flatten weights at energies above the maximum entropy state,
  // and subtract off minimum weight so that our weights don't get out of hand
  void flush_weight_array();

  /*** HISTOGRAM METHODS ***/

  // set canonical weights below some given energy
  void initialize_canonical(double T, int reference=0);

  void initialize_wltmmc(double wl_fmod,
                         double wl_threshold, double wl_cutoff);
  void initialize_satmmc();
  void initialize_samc(bool am_sad = false);
  void initialize_wang_landau(double wl_fmod,
                              double wl_threshold, double wl_cutoff,
                              bool fixed_energy_range);

  void initialize_optimized_ensemble(int first_update_iterations, int oe_update_factor);

  void initialize_simple_flat(int flat_update_factor);

  void initialize_tmi(int version=1);
  void initialize_toe(int version=1);
  void initialize_transitions();

  void initialize_transitions_file(const char *transitions_input_filename);
  void write_transitions_file();
  void write_header(FILE *f);

  double fractional_dos_precision;
  void update_weights_using_transitions(int version, bool energy_range_fixed = false);
  void calculate_weights_using_wltmmc(double wl_fmod,
                                      double wl_threshold, double wl_cutoff,
                                      bool verbose); // added by JP in 2017 for wltmmc.
  void optimize_weights_using_transitions(int version);

  // return fractional error in sample count
  double fractional_sample_error(double T, bool optimistic_sampling);

  double* compute_ln_dos(dos_types dos_type);
  double *compute_walker_density_using_transitions(double *sample_rate = 0);

  int set_min_important_energy(double *ln_dos = 0);
  void set_max_entropy_energy();

  // check whether we are done initializing
  bool finished_initializing(bool be_verbose = false);
  bool reached_iteration_cap();
  // For the following bit, convergence is defined in terms of having
  // 10 pessimistic samples.
  double converged_to_temperature(double *ln_dos) const;
  int converged_to_state() const;

  double estimate_trip_time(int E1, int E2);

  // check whether we may print, to prevent dumping obscene amounts of text into the console
  bool printing_allowed();

  // manual minimum important energies for Wang-Landau
  int default_min_e(){
    if(min_T == 0.2){
      if(N == 5) return 10;
      if(N == 6) return 15;
      if(N == 7) return 20;
      if(N == 8) return 26;
      if(N == 9) return 28;
      if(N == 10) return 37;
      if(N == 11) return 42;
      if(N == 12) return 45;
      if(N == 13) return 58;
      if(N == 14) return 65;
      if(N == 15) return 72;
      if(N == 16) return 76;
      if(N == 17) return 80;
      if(N == 18) return 83;
      if(N == 19) return 88;
      if(N == 20) return 95;
      if(N == 21) return 102;
      if(N == 22) return 108;
      if(N == 23) return 110;
      if(N == 24) return 116;
      if(N == 25) return 123;
      if(N == 26) return 127;
      if(N == 27) return 132;
      if(N == 28) return 140;
      if(N == 29) return 148;
      if(N == 30) return 155;
    }
    printf("\nWe do not know the minimum important energy for the given parameters:");
    printf("N: %d,  min_T: %g\n",N,min_T);
    exit(179);
    return 0;
  }

  sw_simulation(){
    energies_found = 0; // we have not yet found any!
    // seconds per iteration (will be adjusted from actual timing)
    estimated_time_per_iteration = 0.1;
    transitions_filename = 0; // default to NULL pointer here for safety.
    transitions_movie_filename_format = 0; // default to NULL pointer here for safety.
    dos_movie_filename_format = 0; // default to NULL pointer here for safety.
    lnw_movie_filename_format = 0; // default to NULL pointer here for safety.
    transitions_movie_count = 0;
    dos_movie_count = 0;
    lnw_movie_count = 0;
    max_time = 0;
    start_time = clock()/double(CLOCKS_PER_SEC);
    wl_factor = 0.0; // default to no WL method
    sa_t0 = 0.0; // default to no SA method either
    sa_prefactor = 1.0; // default to standard SAMC when using SAMC.
    use_satmmc = false; // default to not using SATMMC.
    use_sad = false; // default to not using SAD MC (think happy thoughts!)
    sa_weight = 0.0;
    use_tmmc = false; // default to not using TMMC for accepting moves.
  };
};

// Return the vector pointing from a to b, accounting for periodic boundaries
vector3d periodic_diff(const vector3d &a, const vector3d  &b,
                          const double len[3], int walls);

// Create and initialize the neighbor tables for all balls (p).
// Returns the maximum number of neighbors that any ball has,
// or -1 if that number is larger than max_neighbors.
int initialize_neighbor_tables(ball *p, int N, double neighborR, int max_neighbors,
                               const double len[3], int walls);

// Find's the neighbors of a by comparing a's position to the center of
// everyone else's neighborsphere, where id is the index of a in p.
void update_neighbors(ball &a, int id, const ball *p, int N,
                      double neighborR, const double len[3], int walls);

// Add ball new_n to the neighbor table of ball id
void add_neighbor(int new_n, ball *p, int id);

// Remove ball old_n from the neighbor table of ball id
void remove_neighbor(int old_n, ball *p, int id);

// Removes p from the neighbor table of anyone neighboring old_p.
// Adds p to the neighbor table of anyone neighboring new_p.
void inform_neighbors(const ball &new_p, const ball &old_p, ball *p);

// Check whether two balls overlap
bool overlap(const ball &a, const ball &b, const double len[3], int walls);

// Check whether ball a overlaps with any of its neighbors in p.
// If count is true, it will return the total number of overlaps, otherwise constant
// it returns 1 if there is at least one overlap, 0 if there are none.
int overlaps_with_any(const ball &a, const ball *p, const double len[3], int walls);

// Count the number of interactions a given ball has
int count_interactions(int id, ball *p, double interaction_scale,
                       double len[3], int walls, int sticky_wall);

// Count the interactions of all the balls
int count_all_interactions(ball *balls, int N, double interaction_scale,
                           double len[3], int walls, int sticky_wall);

// This function finds the maximum number of balls within a given distance
//   distance should be normalized to (divided by) ball radius
int max_balls_within(double radius);
