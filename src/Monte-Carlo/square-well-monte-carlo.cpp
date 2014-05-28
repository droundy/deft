#include <stdio.h>
#include <time.h>
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <popt.h>
#include <sys/stat.h>
#include "handymath.h"
#include "vector3d.h"
#include "Monte-Carlo/square-well.h"

// ------------------------------------------------------------------------------
// Notes on conventions and definitions used
// ------------------------------------------------------------------------------
//
// All coordinates are cartesion.
//
// The coordinates x, y, z are always floating point numbers that refer to real
// locations within the cell.
//
// The coordinates x_i, y_i, z_i are always integers that refer to the bin number
// of the respective coordinate such that x_i refers to a bin of thickness dx
// centered at x.
//
// The symbols e, e_i, and de are used as general coordinates.
//
// Def: If two objects, a and b, are closer than a.R + b.R + neighbor_R + dn,
// then they are neighbors.
//
// Neighbors are used to drastically reduce the number of collision tests needed.
//
// Def: The neighborsphere of an object, a, is the sphere within which
// everything is a neighbor of a.
// Note that this sphere has a well defined center, but it does not have
// a well defined radius unless all obects are circumscribed by spheres of
// the same radius, but this does not affect anything.


// ------------------------------------------------------------------------------
// Global Constants
// ------------------------------------------------------------------------------

const int x = 0;
const int y = 1;
const int z = 2;

// ------------------------------------------------------------------------------
// Functions
// ------------------------------------------------------------------------------

// States how long it's been since last took call.
static void took(const char *name);

// Saves the locations of all balls to a file.
inline void save_locations(const ball *p, int N, const char *fname,
                           const double len[3], const char *comment="");

// The following functions only do anything if debug is true:

// Prints the location and radius of every ball
// As well as whether any overlap or are outside the cell.
inline void print_all(const ball *p, int N, double len[3]);

// Same as print_all, but only prints information for one ball,
// and does not do overlap tests.
inline void print_one(const ball &a, int id, const ball *p, int N,
                      double len[3], int walls);

// Only print those balls that overlap or are outside the cell
// Also prints those they overlap with
inline void print_bad(const ball *p, int N, double len[3], int walls);

// Checks to make sure that every ball is his neighbor's neighbor.
inline void check_neighbor_symmetry(const ball *p, int N);

int main(int argc, const char *argv[]) {
  took("Starting program");
  // ----------------------------------------------------------------------------
  // Define "Constants" -- set from arguments then unchanged
  // ----------------------------------------------------------------------------

  // NOTE: debug can slow things down VERY much
  int debug = false;
  int test_weights = false;
  int print_weights = false;

  int no_weights = false;
  double fix_kT = 0;
  int flat_histogram = false;
  int gaussian_fit = false;
  int walker_weights = false;
  int wang_landau = false;

  double wl_factor = 1;
  double wl_fmod = 2;
  double wl_threshold = 0.1;
  double wl_cutoff = 1e-8;

  double len[3] = {1, 1, 1};
  int walls = 0;
  int wall_dim = 1;
  unsigned long int seed = 0;

  char *dir = new char[1024];
  sprintf(dir, "papers/square-well-liquid/data");
  char *filename = new char[1024];
  sprintf(filename, "default_filename");
  int N = 1000;
  long iterations = 2500000;
  long initialization_iterations = 500000;
  double acceptance_goal = .4;
  double R = 1;
  double well_width = 1.3;
  double ff = 0.3;
  double neighbor_scale = 2;
  double dr = 0.1;
  double de_density = 0.1;
  double de_g = 0.05;
  double max_rdf_radius = 10;
  int totime = 0;
  // scale is not a quite "constant" -- it is adjusted during the initialization
  //  so that we have a reasonable acceptance rate
  double translation_scale = 0.05;

  poptContext optCon;
  // ----------------------------------------------------------------------------
  // Set values from parameters
  // ----------------------------------------------------------------------------
  poptOption optionsTable[] = {
    {"N", '\0', POPT_ARG_INT, &N, 0, "Number of balls to simulate", "INT"},
    {"ww", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &well_width, 0,
     "Ratio of square well width to ball diameter", "DOUBLE"},
    {"ff", '\0', POPT_ARG_DOUBLE, &ff, 0, "Filling fraction. If specified, the "
     "cell dimensions are adjusted accordingly without changing the shape of "
     "the cell"},
    {"walls", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &walls, 0,
     "Number of walled dimensions (dimension order: x,y,z)", "INT"},
    {"initialize", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT,
     &initialization_iterations, 0,
     "Number of iterations to run for initialization", "INT"},
    {"iterations", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &iterations,
     0, "Number of iterations to run for", "INT"},
    {"de_g", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &de_g, 0,
     "Resolution of distribution functions", "DOUBLE"},
    {"dr", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &dr, 0,
     "Differential radius change used in pressure calculation", "DOUBLE"},
    {"de_density", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &de_density, 0, "Resolution of density file", "DOUBLE"},
    {"max_rdf_radius", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &max_rdf_radius, 0, "Set maximum radius for RDF data collection", "DOUBLE"},
    {"lenx", '\0', POPT_ARG_DOUBLE, &len[x], 0,
     "Relative cell size in x dimension", "DOUBLE"},
    {"leny", '\0', POPT_ARG_DOUBLE, &len[y], 0,
     "Relative cell size in y dimension", "DOUBLE"},
    {"lenz", '\0', POPT_ARG_DOUBLE, &len[z], 0,
     "Relative cell size in z dimension", "DOUBLE"},
    {"filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     "Base of output file names", "STRING"},
    {"dir", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &dir, 0,
     "Save directory", "dir"},
    {"neighbor_scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &neighbor_scale, 0, "Ratio of neighbor sphere radius to interaction scale "
     "times ball radius. Drastically reduces collision detections","DOUBLE"},
    {"translation_scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &translation_scale, 0, "Standard deviation for translations of balls, "
     "relative to ball radius", "DOUBLE"},
    {"seed", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0,
     "Seed for the random number generator", "INT"},
    {"acceptance_goal", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &acceptance_goal, 0, "Goal to set the acceptance rate", "DOUBLE"},
    {"nw", '\0', POPT_ARG_NONE, &no_weights, 0, "Don't use weighing method "
     "to get better statistics on low entropy states", "BOOLEAN"},
    {"kT", '\0', POPT_ARG_DOUBLE, &fix_kT, 0, "Use a fixed temperature of kT"
     " rather than adjusted weights", "DOUBLE"},
    {"flat", '\0', POPT_ARG_NONE, &flat_histogram, 0,
     "Use a flat histogram method", "BOOLEAN"},
    {"gaussian", '\0', POPT_ARG_NONE, &gaussian_fit, 0,
     "Use gaussian weights for flat histogram", "BOOLEAN"},
    {"walkers", '\0', POPT_ARG_NONE, &walker_weights, 0,
     "Use a walker optimization weight histogram method", "BOOLEAN"},
    {"wang_landau", '\0', POPT_ARG_NONE, &wang_landau, 0,
     "Use Wang-Landau histogram method", "BOOLEAN"},
    {"wl_factor", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &wl_factor,
     0, "Initial value of Wang-Landau factor", "DOUBLE"},
    {"wl_fmod", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &wl_fmod, 0,
     "Wang-Landau factor modifiction parameter", "DOUBLE"},
    {"wl_threshold", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &wl_threshold, 0, "Threhold for normalized standard deviation in "
     "energy histogram at which to adjust Wang-Landau factor", "DOUBLE"},
    {"wl_cutoff", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &wl_cutoff, 0, "Cutoff for Wang-Landau factor", "DOUBLE"},
    {"time", '\0', POPT_ARG_INT, &totime, 0,
     "Timing of display information (seconds)", "INT"},
    {"R", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &R, 0, "Ball radius (for testing purposes; should always be 1)", "DOUBLE"},
    {"test_weights", '\0', POPT_ARG_NONE, &test_weights, 0,
     "Periodically print weight histogram during initialization", "BOOLEAN"},
    {"debug", '\0', POPT_ARG_NONE, &debug, 0, "Debug mode", "BOOLEAN"},
    POPT_AUTOHELP
    POPT_TABLEEND
  };
  optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTION...]\nNumber of balls and filling "
                         "fraction or cell dimensions are required arguments.");

  int c = 0;
  // go through arguments, set them based on optionsTable
  while((c = poptGetNextOpt(optCon)) >= 0);
  if (c < -1) {
    fprintf(stderr, "\n%s: %s\n", poptBadOption(optCon, 0), poptStrerror(c));
    return 1;
  }
  poptFreeContext(optCon);

  // ----------------------------------------------------------------------------
  // Verify we have reasonable arguments and set secondary parameters
  // ----------------------------------------------------------------------------

  // check that only one method is used
  if(bool(no_weights) + bool(flat_histogram) + bool(gaussian_fit)
     + bool(wang_landau) + bool(walker_weights) + (fix_kT != 0) != 1){
    printf("Exactly one histigram method must be selected!");
    return 254;
  }

  if(walls >= 2){
    printf("Code cannot currently handle walls in more than one dimension.\n");
    return 254;
  }
  if(walls > 3){
    printf("You cannot have walls in more than three dimensions.\n");
    return 254;
  }
  if(well_width < 1){
    printf("Interaction scale should be greater than (or equal to) 1.\n");
    return 254;
  }

  // Adjust cell dimensions for desired filling fraction
  const double fac = R*pow(4.0/3.0*M_PI*N/(ff*len[x]*len[y]*len[z]), 1.0/3.0);
  for(int i = 0; i < 3; i++) len[i] *= fac;
  printf("\nSetting cell dimensions to (%g, %g, %g).\n",
         len[x], len[y], len[z]);
  if (N <= 0 || initialization_iterations < 0 || iterations < 0 || R <= 0 ||
      neighbor_scale <= 0 || dr <= 0 || translation_scale < 0 ||
      len[x] < 0 || len[y] < 0 || len[z] < 0) {
    fprintf(stderr, "\nAll parameters must be positive.\n");
    return 1;
  }
  dr *= R;

  const double eta = (double)N*4.0/3.0*M_PI*R*R*R/(len[x]*len[y]*len[z]);
  if (eta > 1) {
    fprintf(stderr, "\nYou're trying to cram too many balls into the cell. "
            "They will never fit. Filling fraction: %g\n", eta);
    return 7;
  }

  // If a filename was not selected, make a default
  if (strcmp(filename, "default_filename") == 0) {
    char *name_suffix = new char[10];
    char *wall_tag = new char[10];
    if(walls == 0) sprintf(wall_tag,"periodic");
    else if(walls == 1) sprintf(wall_tag,"wall");
    else if(walls == 2) sprintf(wall_tag,"tube");
    else if(walls == 3) sprintf(wall_tag,"box");
    if (fix_kT) {
      sprintf(name_suffix, "-kT%g", fix_kT);
    } else if (no_weights) {
      sprintf(name_suffix, "-nw");
    } else if (flat_histogram) {
      sprintf(name_suffix, "-flat");
    } else if (gaussian_fit) {
      sprintf(name_suffix, "-gaussian");
    } else if (wang_landau) {
      sprintf(name_suffix, "-wang_landau");
    } else if (walker_weights) {
      sprintf(name_suffix, "-walkers");
    } else {
      name_suffix[0] = 0; // set name_suffix to the empty string
    }
    sprintf(filename, "%s-ww%04.2f-ff%04.2f-N%i%s",
            wall_tag, well_width, eta, N, name_suffix);
    printf("\nUsing default file name: ");
    delete[] name_suffix;
    delete[] wall_tag;
  }
  else
    printf("\nUsing given file name: ");
  printf("%s\n",filename);

  printf("------------------------------------------------------------------\n");
  printf("Running %s with parameters:\n", argv[0]);
  for(int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') printf("\n");
    printf("%s ", argv[i]);
  }
  printf("\n");
  if (totime > 0) printf("Timing information will be displayed.\n");
  if (debug) printf("DEBUG MODE IS ENABLED!\n");
  else printf("Debug mode disabled\n");
  printf("------------------------------------------------------------------\n\n");

  // ----------------------------------------------------------------------------
  // Define variables
  // ----------------------------------------------------------------------------

  // translation distance should scale with ball radius
  double translation_distance = translation_scale*R;

  // neighbor radius should scale with radius and interaction scale
  double neighbor_R = neighbor_scale*R*well_width;

  // Find the upper limit to the maximum number of neighbors a ball could have
  int max_neighbors = max_balls_within(2+neighbor_scale*well_width);

  // Energy histogram
  const double interaction_distance = 2*R*well_width;
  const int energy_levels = N/2*max_balls_within(interaction_distance);
  long *energy_histogram = new long[energy_levels]();

  // Walkers
  bool current_walker_plus = false;
  int walker_plus_threshold = 0, walker_minus_threshold = 0;
  long *walkers_plus = new long[energy_levels]();
  long *walkers_total = new long[energy_levels]();

  // Energy weights, state density
  int weight_updates = 0;
  double *ln_energy_weights = new double[energy_levels]();
  if (fix_kT) {
    for(int i = 0; i < energy_levels; i++){
      ln_energy_weights[i] = i/fix_kT;
    }
  }

  // Radial distribution function (RDF) histogram
  long *g_energy_histogram = new long[energy_levels]();
  const int g_bins = round(min(min(min(len[y],len[z]),len[x]),max_rdf_radius)
                           / de_g / 2);
  long **g_histogram = new long*[energy_levels];
  for(int i = 0; i < energy_levels; i++)
    g_histogram[i] = new long[g_bins]();

  // Density histogram
  const int density_bins = round(len[wall_dim]/de_density);
  const double bin_volume = len[x]*len[y]*len[z]/len[wall_dim]*de_density;
  long **density_histogram = new long*[energy_levels];
  for(int i = 0; i < energy_levels; i++)
    density_histogram[i] = new long[density_bins]();

  printf("memory use estimate = %.2g G\n\n",
         8*double((6 + g_bins + density_bins)*energy_levels)/1024/1024/1024);

  ball *balls = new ball[N];
  move_info moves;

  if(totime < 0) totime = 10*N;

  // a guess for the number of iterations to run for initializing the histogram
  const int first_weight_update = energy_levels;

  // hokey, very probably temporary way of running wang_landau
  if(wang_landau) initialization_iterations = N*first_weight_update;

  // Initialize the random number generator with our seed
  random::seed(seed);

  // ----------------------------------------------------------------------------
  // Set up the initial grid of balls
  // ----------------------------------------------------------------------------

  for(int i = 0; i < N; i++) // initialize ball radii
    balls[i].R = R;

  // Balls will be initially placed on a face centered cubic (fcc) grid
  // Note that the unit cells need not be actually "cubic", but the fcc grid will
  //   be stretched to cell dimensions
  const int spots_per_cell = 4; // spots in each fcc periodic unit cell
  const int cells_floor = ceil(N/spots_per_cell); // minimum number of cells
  int cells[3]; // array to contain number of cells in x, y, and z dimensions
  for(int i = 0; i < 3; i++){
    cells[i] = ceil(pow(cells_floor*len[i]*len[i]
                        /(len[(i+1)%3]*len[(i+2)%3]),1.0/3.0));
  }

  // It is usefull to know our cell dimensions
  double cell_width[3];
  for(int i = 0; i < 3; i++) cell_width[i] = len[i]/cells[i];

  // Increase number of cells until all balls can be accomodated
  int total_spots = spots_per_cell*cells[x]*cells[y]*cells[z];
  int i = 0;
  while(total_spots < N) {
    if(cell_width[i%3] <= cell_width[(i+1)%3] &&
       cell_width[(i+1)%3] <= cell_width[(i+2)%3]) {
      cells[i%3] += 1;
      cell_width[i%3] = len[i%3]/cells[i%3];
      total_spots += spots_per_cell*cells[(i+1)%3]*cells[(i+2)%3];
    }
    i++;
  }

  // Define ball positions relative to cell position
  vector3d* offset = new vector3d[4]();
  offset[x] = vector3d(0,cell_width[y],cell_width[z])/2;
  offset[y] = vector3d(cell_width[x],0,cell_width[z])/2;
  offset[z] = vector3d(cell_width[x],cell_width[y],0)/2;

  // Reserve some spots at random to be vacant
  bool *spot_reserved = new bool[total_spots]();
  int p; // Index of reserved spot
  for(int i = 0; i < total_spots-N; i++) {
    p = floor(random::ran()*total_spots); // Pick a random spot index
    if(spot_reserved[p] == false) // If it's not already reserved, reserve it
      spot_reserved[p] = true;
    else // Otherwise redo this index (look for a new spot)
      i--;
  }

  // Place all balls in remaining spots
  int b = 0;
  vector3d cell_pos;
  for(int i = 0; i < cells[x]; i++) {
    for(int j = 0; j < cells[y]; j++) {
      for(int k = 0; k < cells[z]; k++) {
        for(int l = 0; l < 4; l++) {
          if(!spot_reserved[i*(4*cells[z]*cells[y])+j*(4*cells[z])+k*4+l]) {
            balls[b].pos = vector3d(i*cell_width[x],j*cell_width[y],
                                    k*cell_width[z]) + offset[l];
            b++;
          }
        }
      }
    }
  }
  delete[] offset;
  delete[] spot_reserved;
  took("Placement");

  // ----------------------------------------------------------------------------
  // Print info about the initial configuration for troubleshooting
  // ----------------------------------------------------------------------------

  int most_neighbors =
    initialize_neighbor_tables(balls, N, neighbor_R + 2*dr, max_neighbors, len,
                               walls);
  if (most_neighbors < 0) {
    fprintf(stderr, "The guess of %i max neighbors was too low. Exiting.\n",
            max_neighbors);
    return 1;
  }
  printf("Neighbor tables initialized.\n");
  printf("The most neighbors is %i, whereas the max allowed is %i.\n",
         most_neighbors, max_neighbors);

  // ----------------------------------------------------------------------------
  // Make sure initial placement is valid
  // ----------------------------------------------------------------------------

  bool error = false, error_cell = false;
  for(int i = 0; i < N; i++) {
    if (!in_cell(balls[i], len, walls, dr)) {
      error_cell = true;
      error = true;
    }
    for(int j = 0; j < i; j++) {
      if (overlap(balls[i], balls[j], len, walls)) {
        error = true;
        break;
      }
    }
    if (error) break;
  }
  if (error){
    print_bad(balls, N, len, walls);
    printf("Error in initial placement: ");
    if(error_cell) printf("balls placed outside of cell.\n");
    else printf("balls are overlapping.\n");
    return 253;
  }

  fflush(stdout);

  // ----------------------------------------------------------------------------
  // Initialization of cell
  // ----------------------------------------------------------------------------

  double avg_neighbors = 0;
  int interactions =
    count_all_interactions(balls, N, interaction_distance, len, walls);
  double dscale = .1;

  for(long iteration = 1;
      iteration <= initialization_iterations + first_weight_update; iteration++) {
    // ---------------------------------------------------------------
    // Move each ball once
    // ---------------------------------------------------------------
    for(int i = 0; i < N; i++) {
      move_one_ball(i, balls, N, len, walls, neighbor_R, translation_distance,
                    interaction_distance, max_neighbors, dr, &moves,
                    interactions, ln_energy_weights);
      interactions += moves.new_count - moves.old_count;
      energy_histogram[interactions]++;

      if(wang_landau && iteration > first_weight_update){
        ln_energy_weights[interactions] -= wl_factor
          - log(1-1/(energy_histogram[interactions]+1e-10));
        // the first term is from wang-langau
        // the second term is to correct for the fact that wang-landau uses
        //   a ratio of DoS rather than a ratio of weights to calculate the
        //   probability of accepting a move
        // the 1e-10 factor is added for the case energy_histogram[interactions] == 1;
        //   its effect should otherwise be negligibly small
      }
      if(walker_weights){
        walkers_total[interactions]++;
        if(interactions >= walker_minus_threshold)
          current_walker_plus = false;
        else if(interactions <= walker_plus_threshold)
          current_walker_plus = true;
        if(current_walker_plus)
          walkers_plus[interactions]++;
      }
    }
    assert(interactions ==
           count_all_interactions(balls, N, interaction_distance, len, walls));
    // ---------------------------------------------------------------
    // Fine-tune translation scale to reach acceptance goal
    // ---------------------------------------------------------------
    if (iteration % N == 0) {
      const double acceptance_rate =
        (double)(moves.working-moves.working_old)
        /(moves.total-moves.total_old);
      moves.working_old = moves.working;
      moves.total_old = moves.total;
      if (acceptance_rate < acceptance_goal)
        translation_distance /= 1+dscale;
      else
        translation_distance *= 1+dscale;
      // hokey heuristic for tuning dscale
      const double closeness = fabs(acceptance_rate - acceptance_goal)
        / acceptance_rate;
      if(closeness > 0.5) dscale *= 2;
      else if(closeness < dscale*2) dscale/=2;
    }
    // ---------------------------------------------------------------
    // Update weights
    // ---------------------------------------------------------------
    if(!(no_weights || gaussian_fit || fix_kT)){
      if(iteration == first_weight_update){
        flat_hist(energy_histogram, ln_energy_weights, energy_levels);
        weight_updates++;
      } else if((flat_histogram || walker_weights)
                && (iteration > first_weight_update)
                && ((iteration-first_weight_update)
                    % int(first_weight_update*uipow(2,weight_updates)) == 0)){
        printf("Weight update: %d.\n", int(uipow(2,weight_updates)));
        if (flat_histogram)
          flat_hist(energy_histogram, ln_energy_weights, energy_levels);
        else if(walker_weights){
          walker_hist(energy_histogram, ln_energy_weights, energy_levels,
                      walkers_plus, walkers_total, &moves);
        }
        weight_updates++;
        if(test_weights) print_weights = true;
      }
      else if(wang_landau && (iteration == N*first_weight_update)){
        // check whether our histogram is flat enough to update wl_factor
        const double variation = count_variation(energy_histogram, ln_energy_weights, energy_levels);
        if (variation < max(wl_threshold, exp(wl_factor))) {
          wl_factor /= wl_fmod;
          // for wang-landau, only flush energy histogram; keep weights
          flush_arrays(energy_histogram, ln_energy_weights, energy_levels, true);
        }
        // repeat until terminal condition is met
        if(wl_factor > wl_cutoff)
          iteration = first_weight_update;

        // print status text for testing purposes
        printf("\nweight update: %i\n",weight_updates);
        printf("  WL factor: %g\n",wl_factor);
        printf("  count variation: %g\n", variation);
        fflush(stdout);

        weight_updates++;
        if((weight_updates % 10 == 0) && test_weights) print_weights = true;
      }
      // for testing purposes; prints energy histogram and weight array
      if(print_weights){
        char *headerinfo = new char[4096];
        sprintf(headerinfo,
                "# cell dimensions: (%5.2f, %5.2f, %5.2f), walls: %i,"
                " de_density: %g, de_g: %g\n# seed: %li, N: %i, R: %f,"
                " well_width: %g, translation_distance: %g\n"
                "# initialization_iterations: %li, neighbor_scale: %g, dr: %g,"
                " energy_levels: %i\n",
                len[0], len[1], len[2], walls, de_density, de_g, seed, N, R,
                well_width, translation_distance, initialization_iterations,
                neighbor_scale, dr, energy_levels);

        char *countinfo = new char[4096];
        sprintf(countinfo,
                "# iteration: %li, working moves: %li, total moves: %li, "
                "acceptance rate: %g\n",
                iteration, moves.working, moves.total,
                double(moves.working)/moves.total);

        const char *testdir = "test";

        char *w_fname = new char[1024];
        char *e_fname = new char[1024];
        mkdir(dir, 0777); // create save directory
        sprintf(w_fname, "%s/%s", dir, testdir);
        mkdir(w_fname, 0777); // create test directory
        sprintf(w_fname, "%s/%s/%s-w%02i.dat",
                dir, testdir, filename, weight_updates);
        sprintf(e_fname, "%s/%s/%s-E%02i.dat",
                dir, testdir, filename, weight_updates);

        FILE *w_out = fopen(w_fname, "w");
        if (!w_out) {
          fprintf(stderr, "Unable to create %s!\n", w_fname);
          exit(1);
        }
        fprintf(w_out, "%s", headerinfo);
        fprintf(w_out, "%s", countinfo);
        fprintf(w_out, "\n# interactions   value\n");
        for(int i = 0; i < energy_levels; i++)
          fprintf(w_out, "%i  %f\n", i, ln_energy_weights[i]);
        fclose(w_out);

        FILE *e_out = fopen((const char *)e_fname, "w");
        fprintf(e_out, "%s", headerinfo);
        fprintf(e_out, "%s", countinfo);
        fprintf(e_out, "\n# interactions   counts\n");
        for(int i = 0; i < energy_levels; i++)
          fprintf(e_out, "%i  %ld\n",i,energy_histogram[i]);
        fclose(e_out);

        delete[] headerinfo;
        delete[] countinfo;
        delete[] w_fname;
        delete[] e_fname;

        print_weights = false;
      }
    }
    // ---------------------------------------------------------------
    // Print out timing information if desired
    // ---------------------------------------------------------------
    if (totime > 0 && iteration % totime == 0) {
      char *iter = new char[1024];
      sprintf(iter, "%i iterations", totime);
      took(iter);
      delete[] iter;
      printf("Iteration %li, acceptance rate of %g, translation_distance: %g.\n",
             iteration, (double)moves.working/moves.total,
             translation_distance);
      printf("We've had %g updates per kilomove and %g informs per kilomoves, "
             "for %g informs per update.\n",
             1000.0*moves.updates/moves.total,
             1000.0*moves.informs/moves.total,
             (double)moves.informs/moves.updates);
      const long checks_without_tables = moves.total*N;
      int total_neighbors = 0;
      for(int i = 0; i < N; i++) {
        total_neighbors += balls[i].num_neighbors;
        most_neighbors = max(balls[i].num_neighbors, most_neighbors);
      }
      avg_neighbors = double(total_neighbors)/N;
      const long checks_with_tables = moves.total*avg_neighbors
        + N*moves.updates;
      printf("We've done about %.3g%% of the distance calculations we would "
             "have done without tables.\n",
             100.0*checks_with_tables/checks_without_tables);
      printf("The max number of neighbors is %i, whereas the most we've seen is "
             "%i.\n", max_neighbors, most_neighbors);
      printf("Neighbor scale is %g and avg. number of neighbors is %g.\n\n",
             neighbor_scale, avg_neighbors);
      fflush(stdout);
    }
  }
  if(gaussian_fit)
    gaussian_hist(energy_histogram,ln_energy_weights,energy_levels);
  took("Initialization");

  // ----------------------------------------------------------------------------
  // Generate info to put in save files
  // ----------------------------------------------------------------------------

  mkdir(dir, 0777); // create save directory

  char *headerinfo = new char[4096];
  sprintf(headerinfo,
          "# cell dimensions: (%5.2f, %5.2f, %5.2f), walls: %i,"
          " de_density: %g, de_g: %g\n# seed: %li, N: %i, R: %f,"
          " well_width: %g, translation_distance: %g\n"
          "# initialization_iterations: %li, neighbor_scale: %g, dr: %g,"
          " energy_levels: %i\n",
          len[0], len[1], len[2], walls, de_density, de_g, seed, N, R,
          well_width, translation_distance, initialization_iterations,
          neighbor_scale, dr, energy_levels);

  char *e_fname = new char[1024];
  sprintf(e_fname, "%s/%s-E.dat", dir, filename);

  char *w_fname = new char[1024];
  sprintf(w_fname, "%s/%s-lnw.dat", dir, filename);

  char *density_fname = new char[1024];
  sprintf(density_fname, "%s/%s-density-%i.dat", dir, filename, N);

  char *g_fname = new char[1024];
  sprintf(g_fname, "%s/%s-g.dat", dir, filename);

  // ----------------------------------------------------------------------------
  // MAIN PROGRAM LOOP
  // ----------------------------------------------------------------------------

  clock_t output_period = CLOCKS_PER_SEC; // start at outputting every minute
  // top out at one hour interval
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*30;
  clock_t last_output = clock(); // when we last output data

  moves.total = 0;
  moves.working = 0;

  // Reset energy histogram
  for(int i = 0; i < energy_levels; i++)
    energy_histogram[i] = 0;

  for(long iteration = 1; iteration <= iterations; iteration++) {
    // ---------------------------------------------------------------
    // Move each ball once, add to energy histogram
    // ---------------------------------------------------------------
    for(int i = 0; i < N; i++) {
      move_one_ball(i, balls, N, len, walls, neighbor_R, translation_distance,
                    interaction_distance, max_neighbors, dr, &moves,
                    interactions, ln_energy_weights);
      interactions += moves.new_count - moves.old_count;
      energy_histogram[interactions]++;
    }
    assert(interactions ==
           count_all_interactions(balls, N, interaction_distance, len, walls));
    // ---------------------------------------------------------------
    // Add data to density and RDF histograms
    // ---------------------------------------------------------------
    // Density histogram
    if(walls){
      for(int i = 0; i < N; i++){
        density_histogram[interactions]
          [int(floor(balls[i].pos[wall_dim]/de_density))] ++;
      }
    }

    // RDF
    if(!walls){
      g_energy_histogram[interactions]++;
      for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
          if(i != j){
            const vector3d r = periodic_diff(balls[i].pos, balls[j].pos, len,
                                             walls);
            const int r_i = floor(r.norm()/de_g);
            if(r_i < g_bins) g_histogram[interactions][r_i]++;
          }
        }
      }
    }
    // ---------------------------------------------------------------
    // Save to file
    // ---------------------------------------------------------------
    const clock_t now = clock();
    if ((now - last_output > output_period) || iteration == iterations) {
      last_output = now;
      assert(last_output);
      if (output_period < max_output_period/2) output_period *= 2;
      else if (output_period < max_output_period)
        output_period = max_output_period;
      const double secs_done = double(now)/CLOCKS_PER_SEC;
      const int seconds = int(secs_done) % 60;
      const int minutes = int(secs_done / 60) % 60;
      const int hours = int(secs_done / 3600) % 24;
      const int days = int(secs_done / 86400);
      printf("Saving data after %i days, %02i:%02i:%02i, %li iterations "
             "complete.\n", days, hours, minutes, seconds, iteration);
      fflush(stdout);

      char *countinfo = new char[4096];
      sprintf(countinfo,
              "# iteration: %li, working moves: %li, total moves: %li, "
              "acceptance rate: %g\n",
              iteration, moves.working, moves.total,
              double(moves.working)/moves.total);

      // Save energy histogram
      FILE *e_out = fopen((const char *)e_fname, "w");
      fprintf(e_out, "%s", headerinfo);
      fprintf(e_out, "%s", countinfo);
      fprintf(e_out, "\n# interactions   counts\n");
      for(int i = 0; i < energy_levels; i++)
        if(energy_histogram[i] != 0)
          fprintf(e_out, "%i  %ld\n",i,energy_histogram[i]);
      fclose(e_out);

      // Save weights histogram
      FILE *w_out = fopen((const char *)w_fname, "w");
      fprintf(w_out, "%s", headerinfo);
      fprintf(w_out, "%s", countinfo);
      fprintf(w_out, "\n# interactions   ln(weight)\n");
      for(int i = 0; i < energy_levels; i++)
        if(energy_histogram[i] != 0){
          fprintf(w_out, "%i  %g\n",i,ln_energy_weights[i]
                  + (wang_landau ? log(energy_histogram[i]) : 0));
        }
      fclose(w_out);

      // Save RDF
      if(!walls){
        FILE *g_out = fopen((const char *)g_fname, "w");
        fprintf(g_out, "%s", headerinfo);
        fprintf(g_out, "%s", countinfo);
        fprintf(g_out, "\n# data table containing values of g");
        fprintf(g_out, "\n# first column reserved for specifying energy level");
        fprintf(g_out, "\n# column number rn (starting from the second column, "
                "counting from zero) corresponds to radius r given by "
                "r = (rn + 0.5) * de_g");
        const double density = N/len[x]/len[y]/len[z];
        const double total_vol = len[x]*len[y]*len[z];
        for(int i = 0; i < energy_levels; i++){
          if(g_histogram[i][g_bins-1] > 0){ // if we have RDF data at this energy
            fprintf(g_out, "\n%i",i);
            for(int r_i = 0; r_i < g_bins; r_i++) {
              const double probability = (double)g_histogram[i][r_i]
                / g_energy_histogram[i];
              const double r = (r_i + 0.5) * de_g;
              const double shell_vol =
                4.0/3.0*M_PI*(uipow(r+de_g/2, 3) - uipow(r-de_g/2, 3));
              const double n2 = probability/total_vol/shell_vol;
              const double g = n2/sqr(density);
              fprintf(g_out, " %8.5f", g);
            }
          }
        }
        fclose(g_out);
      }

      // Saving density data
      if(walls){
        FILE *densityout = fopen((const char *)density_fname, "w");
        fprintf(densityout, "%s", headerinfo);
        fprintf(densityout, "%s", countinfo);
        fprintf(densityout, "\n# data table containing densities in slabs "
                "(bins) of thickness de_density away from a wall");
        fprintf(densityout, "\n# row number corresponds to energy level");
        fprintf(densityout, "\n# column number dn (counting from zero) "
                "corresponds to distance d from wall given by "
                "d = (dn + 0.5) * de_density");
        for(int i = 0; i < energy_levels; i++){
          fprintf(densityout, "\n");
          for(int r_i = 0; r_i < density_bins; r_i++) {
            const double bin_density =
              (double)density_histogram[i][r_i]
              *N/energy_histogram[i]/bin_volume;
            fprintf(densityout, "%8.5f ", bin_density);
          }
        }
        fclose(densityout);
      }

      delete[] countinfo;
    }
  }
  // ----------------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // ----------------------------------------------------------------------------

  delete[] balls;
  delete[] g_histogram;
  delete[] density_histogram;

  delete[] headerinfo;
  return 0;
}
// ------------------------------------------------------------------------------
// END OF MAIN
// ------------------------------------------------------------------------------

inline void print_all(const ball *p, int N) {
  for (int i = 0; i < N; i++) {
    char *pos = new char[1024];
    p[i].pos.tostr(pos);
    printf("%4i: R: %4.2f, %i neighbors: ", i, p[i].R, p[i].num_neighbors);
    for(int j = 0; j < min(10, p[i].num_neighbors); j++)
      printf("%i ", p[i].neighbors[j]);
    if (p[i].num_neighbors > 10)
      printf("...");
    printf("\n      pos:          %s\n", pos);
    delete[] pos;
  }
  printf("\n");
  fflush(stdout);
}

inline void print_one(const ball &a, int id, const ball *p, int N,
                      double len[3], int walls) {
  char *pos = new char[1024];
  a.pos.tostr(pos);
  printf("%4i: R: %4.2f, %i neighbors: ", id, a.R, a.num_neighbors);
  for(int j=0; j<min(10, a.num_neighbors); j++)
    printf("%i ", a.neighbors[j]);
  if (a.num_neighbors > 10)
    printf("...");
  printf("\n      pos:          %s\n", pos);
  for (int j=0; j<N; j++) {
    if (j != id && overlap(a, p[j], len, walls)) {
      p[j].pos.tostr(pos);
      printf("\t  Overlaps with %i", j);
      printf(": %s\n", pos);
    }
  }
  delete[] pos;
  printf("\n");
  fflush(stdout);
}

inline void print_bad(const ball *p, int N, double len[3], int walls) {
  for (int i = 0; i < N; i++) {
    bool incell = in_cell(p[i], len, walls);
    bool overlaps = false;
    for (int j = 0; j < i; j++) {
      if (overlap(p[i], p[j], len, walls)) {
        overlaps = true;
        break;
      }
    }
    if (!incell || overlaps) {
      char *pos = new char[1024];
      p[i].pos.tostr(pos);
      printf("%4i: %s R: %4.2f\n", i, pos, p[i].R);
      if (!incell)
        printf("\t  Outside cell!\n");
      for (int j = 0; j < i; j++) {
        if (overlap(p[i], p[j], len, walls)) {
          p[j].pos.tostr(pos);
          printf("\t  Overlaps with %i", j);
          printf(": %s\n", pos);
        }
      }
      delete[] pos;
    }
  }
  fflush(stdout);
}

inline void check_neighbor_symmetry(const ball *p, int N) {
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < p[i].num_neighbors; j++) {
      const int k = p[i].neighbors[j];
      bool is_neighbor = false;
      for (int l = 0; l < p[k].num_neighbors; l++) {
        if (p[k].neighbors[l] == i) {
          is_neighbor = true;
          break;
        }
      }
      if(!is_neighbor) {
        printf("NEIGHBOR TABLE ERROR: %i has %i as a neighbor, but %i does "
               "not reciprocate!!!\n", i, k, k);
      }
    }
  }
}

static void took(const char *name) {
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
}

void save_locations(const ball *p, int N, const char *fname, const double len[3],
                    const char *comment) {
  FILE *out = fopen((const char *)fname, "w");
  fprintf(out, "# %s\n", comment);
  fprintf(out, "%g %g %g\n", len[x], len[y], len[z]);
  for(int i = 0; i < N; i++) {
    fprintf(out, "%6.2f %6.2f %6.2f ", p[i].pos[x], p[i].pos[y], p[i].pos[z]);
    fprintf(out, "\n");
  }
  fclose(out);
}
