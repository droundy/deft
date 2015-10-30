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

// Tests validity of a shrunken version of the cell
static bool overlap_in_small_cell(sw_simulation &sw, double scaling_factor);

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

// Checks to make sure that every ball is his neighbor's neighbor.
inline void check_neighbor_symmetry(const ball *p, int N);

// Check whether it is a reasonable time to save data according to our
// heuristic.
bool time_to_save();

int main(int argc, const char *argv[]) {
  took("Starting program");
  // ----------------------------------------------------------------------------
  // Define "Constants" -- set from arguments then unchanged
  // ----------------------------------------------------------------------------

  // NOTE: debug can slow things down VERY much
  int debug = false;

  sw_simulation sw;

  sw.sticky_wall = 0;
  sw.len[0] = sw.len[1] = sw.len[2] = 1;
  sw.walls = 0;
  sw.N = 10;
  sw.translation_scale = 0.05;

  unsigned long int seed = 0;

  char *data_dir = new char[1024];
  sprintf(data_dir, "free-energy-data");
  char *filename = new char[1024];
  sprintf(filename, "default_filename");
  char *filename_suffix = new char[1024];
  sprintf(filename_suffix, "default_filename_suffix");
  long simulation_iterations = 1000000;
  double acceptance_goal = .4;
  double R = 1;
  const double well_width = 1;
  double ff_small = -1;
  double neighbor_scale = 2;
  double de_g = 0.05;
  double max_rdf_radius = 10;
  int totime = 0;

  poptContext optCon;
  // ----------------------------------------------------------------------------
  // Set values from parameters
  // ----------------------------------------------------------------------------
  poptOption optionsTable[] = {
    {"N", '\0', POPT_ARG_INT, &sw.N, 0, "Number of balls to simulate", "INT"},
    {"ff_small", '\0', POPT_ARG_DOUBLE, &ff_small, 0, "Small filling fraction. If specified, "
     "This sets the desired filling fraction of the shrunk cell. Otherwise it defaults to ff."},
    {"walls", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &sw.walls, 0,
     "Number of walled dimensions (dimension order: x,y,z)", "INT"},
    {"iterations", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &simulation_iterations,
     0, "Number of iterations for which to run the simulation", "INT"},
    {"de_g", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &de_g, 0,
     "Resolution of distribution functions", "DOUBLE"},
    {"max_rdf_radius", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &max_rdf_radius, 0, "Set maximum radius for RDF data collection", "DOUBLE"},
    {"lenx", '\0', POPT_ARG_DOUBLE, &sw.len[x], 0,
     "Relative cell size in x dimension", "DOUBLE"},
    {"leny", '\0', POPT_ARG_DOUBLE, &sw.len[y], 0,
     "Relative cell size in y dimension", "DOUBLE"},
    {"lenz", '\0', POPT_ARG_DOUBLE, &sw.len[z], 0,
     "Relative cell size in z dimension", "DOUBLE"},
    {"filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     "Base of output file names", "STRING"},
    {"filename_suffix", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
     &filename_suffix, 0, "Output file name suffix", "STRING"},
    {"data_dir", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &data_dir, 0,
     "Directory in which to save data", "data_dir"},
    {"neighbor_scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &neighbor_scale, 0, "Ratio of neighbor sphere radius to interaction scale "
     "times ball radius. Drastically reduces collision detections","DOUBLE"},
    {"translation_scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &sw.translation_scale, 0, "Standard deviation for translations of balls, "
     "relative to ball radius", "DOUBLE"},
    {"seed", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0,
     "Seed for the random number generator", "INT"},
    {"acceptance_goal", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &acceptance_goal, 0, "Goal to set the acceptance rate", "DOUBLE"},
    {"time", '\0', POPT_ARG_INT, &totime, 0,
     "Timing of display information (seconds)", "INT"},
    {"R", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &R, 0, "Ball radius (for testing purposes; should always be 1)", "DOUBLE"},
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

  if(sw.walls >= 2){
    printf("Code cannot currently handle walls in more than one dimension.\n");
    return 254;
  }
  if(sw.walls > 3){
    printf("You cannot have walls in more than three dimensions.\n");
    return 254;
  }
  if(well_width < 1){
    printf("Interaction scale should be greater than (or equal to) 1.\n");
    return 254;
  }

  if (ff_small == -1) {
    printf("You must specify --ff_small!\n    Silly you!\n");
    exit(1);
  }

  if (ff_small != 0) {
    // The user specified a filling fraction, so we must make it so!
    const double volume = 4*M_PI/3*R*R*R*sw.N/ff_small;
    const double min_cell_width = 2*sqrt(2)*R; // minimum cell width
    const int numcells = (sw.N+3)/4; // number of unit cells we need
    const int max_cubic_width
      = pow(volume/min_cell_width/min_cell_width/min_cell_width, 1.0/3);
    if (max_cubic_width*max_cubic_width*max_cubic_width >= numcells) {
      // We can get away with a cubic cell, so let's do so.  Cubic
      // cells are nice and comfortable!
      sw.len[x] = sw.len[y] = sw.len[z] = pow(volume, 1.0/3);
    } else {
      // A cubic cell won't work with our initialization routine, so
      // let's go with a lopsided cell that should give us something
      // that will work.
      int xcells = int( pow(numcells, 1.0/3) );
      int cellsleft = (numcells + xcells - 1)/xcells;
      int ycells = int( sqrt(cellsleft) );
      int zcells = (cellsleft + ycells - 1)/ycells;

      // The above should give a zcells that is largest, followed by
      // ycells and then xcells.  Thus we make the lenz just small
      // enough to fit these cells, and so on, to make the cell as
      // close to cubic as possible.
      sw.len[z] = zcells*min_cell_width;
      if (xcells == ycells) {
        sw.len[x] = sw.len[y] = sqrt(volume/sw.len[z]);
      } else {
        sw.len[y] = min_cell_width*ycells;
        sw.len[x] = volume/sw.len[y]/sw.len[z];
      }
      printf("Using lopsided %d x %d x %d cell (total goal %d)\n",
             xcells, ycells, zcells, numcells);
    }
  }


  printf("\nSetting cell dimensions to (%g, %g, %g).\n",
         sw.len[x], sw.len[y], sw.len[z]);
  printf("\nFilling fraction of small cell is %g", ff_small);
  if (sw.N <= 0 || simulation_iterations < 0 || R <= 0 ||
      neighbor_scale <= 0 || sw.translation_scale < 0 ||
      sw.len[x] < 0 || sw.len[y] < 0 || sw.len[z] < 0) {
    fprintf(stderr, "\nAll parameters must be positive.\n");
    return 1;
  }

  const double eta = (double)sw.N*4.0/3.0*M_PI*R*R*R/(sw.len[x]*sw.len[y]*sw.len[z]);
  if (eta > 1) {
    fprintf(stderr, "\nYou're trying to cram too many balls into the cell. "
            "They will never fit. Filling fraction (small): %g\n", eta);
    return 7;
  }

  // If a filename was not selected, make a default
  if (strcmp(filename, "default_filename") == 0) {
    char *wall_tag = new char[100];
    if(sw.walls == 0) sprintf(wall_tag,"periodic");
    else if(sw.walls == 1) sprintf(wall_tag,"wall");
    else if(sw.walls == 2) sprintf(wall_tag,"tube");
    else if(sw.walls == 3) sprintf(wall_tag,"box");
    sprintf(filename, "%s-ww%04.2f-ff_small%04.2f-N%i",
            wall_tag, well_width, eta, sw.N);
    printf("\nUsing default file name: ");
    delete[] wall_tag;
  }
  else
    printf("\nUsing given file name: ");
  // If a filename suffix was specified, add it
  if (strcmp(filename_suffix, "default_filename_suffix") != 0)
    sprintf(filename, "%s-%s", filename, filename_suffix);
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
  // Define sw_simulation variables
  // ----------------------------------------------------------------------------

  sw.iteration = 0; // start at zeroeth iteration
  sw.max_entropy_state = 0;
  sw.min_energy_state = 0;
  sw.energy = 0;

  // translation distance should scale with ball radius
  sw.translation_scale *= R;

  // neighbor radius should scale with radius and interaction scale
  sw.neighbor_R = neighbor_scale*R*well_width;

  // Find the upper limit to the maximum number of neighbors a ball could have
  sw.max_neighbors = max_balls_within(2+neighbor_scale*well_width);

  // Energy histogram and weights
  sw.interaction_distance = 2*R*well_width;
  sw.energy_levels = 1; // hard spheres can only have one energy!
  sw.energy_histogram = new long[sw.energy_levels]();
  sw.ln_energy_weights = new double[sw.energy_levels]();

  // Observed and sampled energies
  sw.pessimistic_observation = new bool[sw.energy_levels]();
  sw.pessimistic_samples = new long[sw.energy_levels]();
  sw.optimistic_samples = new long[sw.energy_levels]();

  // Transitions from one energy to another
  sw.biggest_energy_transition = max_balls_within(sw.interaction_distance);
  sw.transitions_table =
    new long[sw.energy_levels*(2*sw.biggest_energy_transition+1)]();

  // Walker histograms
  sw.walkers_up = new long[sw.energy_levels]();


  // Radial distribution function (RDF) histogram
  long *g_energy_histogram = new long[sw.energy_levels]();
  const int g_bins = round(min(min(min(sw.len[y],sw.len[z]),sw.len[x]),max_rdf_radius)
                           / de_g / 2);
  long **g_histogram = new long*[sw.energy_levels];
  for(int i = 0; i < sw.energy_levels; i++)
    g_histogram[i] = new long[g_bins]();

  printf("memory use estimate = %.2g G\n\n",
         8*double((6 + g_bins)*sw.energy_levels)/1024/1024/1024);

  sw.balls = new ball[sw.N];

  if(totime < 0) totime = 10*sw.N;

  // Initialize the random number generator with our seed
  random::seed(seed);

  // ----------------------------------------------------------------------------
  // Set up the initial grid of balls
  // ----------------------------------------------------------------------------

  for(int i = 0; i < sw.N; i++) // initialize ball radii
    sw.balls[i].R = R;

  // Balls will be initially placed on a face centered cubic (fcc) grid
  // Note that the unit cells need not be actually "cubic", but the fcc grid will
  //   be stretched to cell dimensions
  const double min_cell_width = 2*sqrt(2)*R; // minimum cell width
  const int spots_per_cell = 4; // spots in each fcc periodic unit cell
  int cells[3]; // array to contain number of cells in x, y, and z dimensions
  for(int i = 0; i < 3; i++){
    cells[i] = int(sw.len[i]/min_cell_width); // max number of cells that will fit
  }

  // It is usefull to know our cell dimensions
  double cell_width[3];
  for(int i = 0; i < 3; i++) cell_width[i] = sw.len[i]/cells[i];

  // If we made our cells to small, return with error
  for(int i = 0; i < 3; i++){
    if(cell_width[i] < min_cell_width){
      printf("Placement cell size too small: (%g,  %g,  %g) coming from (%g, %g, %g)\n",
             cell_width[0],cell_width[1],cell_width[2],
             sw.len[0], sw.len[1], sw.len[2]);
      printf("Minimum allowed placement cell width: %g\n",min_cell_width);
      printf("Total simulation cell dimensions: (%g,  %g,  %g)\n",
             sw.len[0],sw.len[1],sw.len[2]);
      printf("Fixing the chosen ball number, filling fractoin, and relative\n"
             "  simulation cell dimensions simultaneously does not appear to be possible\n");
      return 176;
    }
  }

  // Define ball positions relative to cell position
  vector3d* offset = new vector3d[4]();
  offset[x] = vector3d(0,cell_width[y],cell_width[z])/2;
  offset[y] = vector3d(cell_width[x],0,cell_width[z])/2;
  offset[z] = vector3d(cell_width[x],cell_width[y],0)/2;

  // Reserve some spots at random to be vacant
  const int total_spots = spots_per_cell*cells[x]*cells[y]*cells[z];
  bool *spot_reserved = new bool[total_spots]();
  int p; // Index of reserved spot
  for(int i = 0; i < total_spots-sw.N; i++) {
    p = floor(random::ran()*total_spots); // Pick a random spot index
    if(spot_reserved[p] == false) // If it's not already reserved, reserve it
      spot_reserved[p] = true;
    else // Otherwise redo this index (look for a new spot)
      i--;
  }

  // Place all balls in remaining spots
  int b = 0;
  for(int i = 0; i < cells[x]; i++) {
    for(int j = 0; j < cells[y]; j++) {
      for(int k = 0; k < cells[z]; k++) {
        for(int l = 0; l < 4; l++) {
          if(!spot_reserved[i*(4*cells[z]*cells[y])+j*(4*cells[z])+k*4+l]) {
            sw.balls[b].pos = vector3d(i*cell_width[x],j*cell_width[y],
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

  {
    int most_neighbors =
      initialize_neighbor_tables(sw.balls, sw.N, sw.neighbor_R,
                                 sw.max_neighbors, sw.len, sw.walls);
    if (most_neighbors < 0) {
      fprintf(stderr, "The guess of %i max neighbors was too low. Exiting.\n",
              sw.max_neighbors);
      return 1;
    }
    printf("Neighbor tables initialized.\n");
    printf("The most neighbors is %i, whereas the max allowed is %i.\n",
           most_neighbors, sw.max_neighbors);
  }


  fflush(stdout);

  // --------------------------------------------------------------------------
  // end initilization routine.
  // --------------------------------------------------------------------------

  // ----------------------------------------------------------------------------
  // Generate info to put in save files
  // ----------------------------------------------------------------------------

  mkdir(data_dir, 0777); // create save directory

  char *headerinfo = new char[4096];
  sprintf(headerinfo,
          "# cell dimensions: (%g, %g, %g)\n"
          "# walls: %i\n"
          "# de_g: %g\n"
          "# seed: %li\n"
          "# N: %i\n"
          "# R: %f\n"
          "# well_width: %g\n"
          "# neighbor_scale: %g\n"
          "# ff_small: %g\n",
          sw.len[0], sw.len[1], sw.len[2], sw.walls, de_g, seed, sw.N, R,
          well_width, neighbor_scale, ff_small);

  char *g_fname = new char[1024];
  sprintf(g_fname, "%s/%s-g.dat", data_dir, filename);

  took("Initialization");

  // ----------------------------------------------------------------------------
  // MAIN PROGRAM LOOP
  // ----------------------------------------------------------------------------

  sw.moves.total = 0;
  sw.moves.working = 0;
  sw.iteration = 0;

  // tracking for free energy
  int total_checks_of_small_cell = 0;
  int total_valid_small_checks = 0;
  int total_failed_small_checks = 0;

  // Reset energy histogram and sample counts
  for(int i = 0; i < sw.energy_levels; i++){
    sw.energy_histogram[i] = 0;
    sw.pessimistic_samples[i] = 0;
    sw.optimistic_samples[i] = 0;
  }

  while(sw.iteration <= simulation_iterations) {
    sw.iteration++;   // not calling move_a_ball(), so end_move_updates() is never called
    // ---------------------------------------------------------------
    // Move each ball once, add to energy histogram
    // ---------------------------------------------------------------
    for(int i = 0; i < sw.N; i++){
      // get a new location for each ball like
      for(int j = 0; j < 3; j++){
        sw.balls[i].pos[j] = sw.len[j] * random::ran();
      }
    }

    total_checks_of_small_cell++;

    double scaling_factor = 1;  // no "scaling" in infinite case
    if(overlap_in_small_cell(sw,  scaling_factor)){
      total_failed_small_checks++;
      if(debug){printf("false\n");}
    }
    else{
      total_valid_small_checks++;
      if(debug){printf("true\n");}
    }

    // ---------------------------------------------------------------
    // Add data to RDF histogram
    // ---------------------------------------------------------------
    if(!sw.walls){
      g_energy_histogram[0]++;
      for(int i = 0; i < sw.N; i++){
        for(int j = 0; j < sw.N; j++){
          if(i != j){
            const vector3d r = periodic_diff(sw.balls[i].pos, sw.balls[j].pos, sw.len,
                                             sw.walls);
            const int r_i = floor(r.norm()/de_g);
            if(r_i < g_bins) g_histogram[0][r_i]++;
          }
        }
      }
    }
    // ---------------------------------------------------------------
    // Save to file
    // ---------------------------------------------------------------

    if (time_to_save() || sw.iteration == simulation_iterations) {
      const clock_t now = clock();
      const double secs_done = double(now)/CLOCKS_PER_SEC;
      const int seconds = int(secs_done) % 60;
      const int minutes = int(secs_done / 60) % 60;
      const int hours = int(secs_done / 3600) % 24;
      const int days = int(secs_done / 86400);
      printf("Saving data after %i days, %02i:%02i:%02i, %li iterations "
             "complete.\n", days, hours, minutes, seconds, sw.iteration);
      fflush(stdout);

      char *countinfo = new char[4096];
      sprintf(countinfo,
              "# iterations: %li\n"
              "# total checks of small cell: %i\n"
              "# total failed small checks: %i\n"
              "# total valid small checks: %i\n\n",
              sw.iteration,
              total_checks_of_small_cell, total_failed_small_checks,
              total_valid_small_checks);

      // Save RDF
      if(!sw.walls){
        FILE *g_out = fopen((const char *)g_fname, "w");
        fprintf(g_out, "%s", headerinfo);
        fprintf(g_out, "%s", countinfo);
        fprintf(g_out, "# data table containing values of g "
                "(i.e. radial distribution function)\n"
                "# first column reserved for specifying energy level\n"
                "# column number r_n (starting from the second column, "
                "counting from zero) corresponds to radius r given by "
                "r = (r_n + 0.5) * de_g\n");
        const double density = sw.N/sw.len[x]/sw.len[y]/sw.len[z];
        const double total_vol = sw.len[x]*sw.len[y]*sw.len[z];
        for(int i = 0; i < sw.energy_levels; i++){
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

      delete[] countinfo;
    }
  }
  // ----------------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // ----------------------------------------------------------------------------

  for (int i=0; i<sw.N; i++) {
    delete[] sw.balls[i].neighbors;
  }
  delete[] sw.balls;
  delete[] sw.ln_energy_weights;
  delete[] sw.energy_histogram;

  delete[] sw.transitions_table;

  delete[] sw.walkers_up;

  delete[] sw.pessimistic_observation;
  delete[] sw.pessimistic_samples;
  delete[] sw.optimistic_samples;

  for (int i = 0; i < sw.energy_levels; i++) {
    delete[] g_histogram[i];
  }
  delete[] g_histogram;
  delete[] g_energy_histogram;

  delete[] headerinfo;
  delete[] g_fname;

  delete[] data_dir;
  delete[] filename;
  delete[] filename_suffix;

  return 0;
}
// ------------------------------------------------------------------------------
// END OF MAIN
// ------------------------------------------------------------------------------

static bool overlap_in_small_cell(sw_simulation &sw, double scaling_factor){
  double scaled_len[3];

  for(int i=0; i < 3; i++){
    scaled_len[i] = scaling_factor * sw.len[i];
  }

  for(int i=0; i<sw.N; i++){
    for (int j=i+1; j<sw.N; j++) {
      // copy pasting from overlap() in square-well.cpp for now
      // contemplated adding a scaling parametor to overlap(), but decided on this for now.
      const vector3d ab = periodic_diff(scaling_factor * sw.balls[i].pos, scaling_factor * sw.balls[j].pos, scaled_len, sw.walls);
      if (ab.normsquared() < sqr(sw.balls[i].R + sw.balls[j].R)) {
        return true;
      }
    }
  }
  return false;
}

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

bool time_to_save() {
  clock_t output_period = CLOCKS_PER_SEC; // start at outputting every second
  // top out at one hour interval
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*60;
  static clock_t last_save_time = 0;

  static int iterations = 0;
  static int how_often = 1;
  // clock can be expensive under fac, so this is a heuristic to
  // reduce our use of it.
  if (++iterations % how_often == 0) {
    const clock_t time_now = clock();
    if(time_now-last_save_time > output_period){

      if (output_period < max_output_period/2) output_period *= 2;
      else if (output_period < max_output_period)
        output_period = max_output_period;

      how_often = 1+ iterations/3; // our simple heuristic
      last_save_time = time_now;
      iterations = 0;
      // flushing occasionally will be no problem and can be helpful
      // if we forget.
      fflush(stdout);
      return true;
    }
  }
  return false;
}
