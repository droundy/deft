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
#include "Monte-Carlo/square-well.h"

#include "version-identifier.h"

// ------------------------------------------------------------------------------
// Notes on conventions and definitions use
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
// Def: If two objects, a and b, are closer than a.R + b.R + neighbor_R,
// then they are neighbors.
//
// Neighbors are used to drastically reduce the number of collision tests needed.
//
// Def: The neighbor sphere of an object, a, is the sphere within which
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
  printf("version: %s\n",version_identifier());
  // ----------------------------------------------------------------------------
  // Define "Constants" -- set from arguments then unchanged
  // ----------------------------------------------------------------------------

  sw_simulation sw;

  // NOTE: debug can slow things down VERY much
  int debug = false;

  int no_weights = false;
  double fix_kT = 0;
  int tmi = false;
  int tmi_version = 1;
  int toe = false;
  int tmmc = false;
  int oetmmc = false;
  int wang_landau = false;
  int vanilla_wang_landau = false;
  int wltmmc = false;
  int simple_flat = false;
  int optimized_ensemble = false;
  int transition_override = false;

  char *transitions_input_filename = new char[1024];
  sprintf(transitions_input_filename, "none");
  int golden = false;

  // Tuning factors
  int oe_update_factor = 2;
  int flat_update_factor = 2;
  sw.min_important_energy = 0;

  /* Do not change these here! They are taken directly from the WL paper.
     If you want to change the WL parameters, run this code with appropriate arguments */
  double wl_factor = 1;
  double wl_fmod = 2;
  double wl_threshold = 0.8;
  double wl_cutoff = 1e-8;

  // end conditions
  int default_pessimistic_min_samples = 10;
  int default_optimistic_min_samples = 200;
  double default_pessimistic_sample_error = 0.2;
  double default_optimistic_sample_error = 0.1;
  double default_flatness = 0.1;
  int optimistic_sampling = true; // popt requires it to be an int, not a bool

  int tmmc_min_samples = 2000;

  // some miscellaneous default or dummy simulation parameters
  sw.len[0] = sw.len[1] = sw.len[2] = 0;
  sw.walls = 0;
  sw.sticky_wall = 0;
  sw.well_width = 1.3;
  sw.filling_fraction = 0.3;
  sw.N = 10;
  sw.translation_scale = 0;
  sw.fractional_dos_precision = 1e-7;
  sw.end_condition = none;
  sw.min_T = 0.2;
  sw.min_samples = 0;
  sw.sample_error = 0;
  sw.flatness = 0;
  sw.init_iters = 0;

  const int wall_dim = 0; // this is hard-coded all over the place:  you can't change it!
  unsigned long int seed = 0;

  char *data_dir = new char[1024];
  sprintf(data_dir,"none");
  char *default_data_dir = new char[1024];
  sprintf(default_data_dir, "papers/histogram/data");
  char *filename = new char[1024];
  sprintf(filename, "none");
  char *filename_suffix = new char[1024];
  sprintf(filename_suffix, "none");
  long simulation_iterations = 3*1e6;
  long simulation_round_trips = 0;
  double acceptance_goal = .4;
  double R = 1;
  double neighbor_scale = 2;
  double de_density = 0.1;
  double de_g = 0.05;
  double max_rdf_radius = 10;
  // scale is not universally constant -- it is adjusted during initialization
  //  so that we have a reasonable acceptance rate

  poptContext optCon;

  // ----------------------------------------------------------------------------
  // Parse input options
  // ----------------------------------------------------------------------------

  poptOption optionsTable[] = {

    /*** FLUID IDENTITY ***/

    {"N", '\0', POPT_ARG_INT, &sw.N, 0, "Number of balls to simulate", "INT"},
    {"ww", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &sw.well_width, 0,
     "Ratio of square well width to ball diameter", "DOUBLE"},
    {"ff", '\0', POPT_ARG_DOUBLE, &sw.filling_fraction, 0, "If specified, the "
     "cell dimensions are adjusted accordingly without changing the shape of the cell"},
    {"walls", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &sw.walls, 0,
     "Number of walled dimensions (dimension order: x,y,z)", "INT"},
    {"sticky-wall", '\0', POPT_ARG_NONE, &sw.sticky_wall, 0, "Make one wall sticky", 0},

    /*** RELATIVE SIZE OF CELL DIMENSIONS ***/

    {"lenx", '\0', POPT_ARG_DOUBLE, &sw.len[x], 0,
     "Relative cell size in x dimension", "DOUBLE"},
    {"leny", '\0', POPT_ARG_DOUBLE, &sw.len[y], 0,
     "Relative cell size in y dimension", "DOUBLE"},
    {"lenz", '\0', POPT_ARG_DOUBLE, &sw.len[z], 0,
     "Relative cell size in z dimension", "DOUBLE"},

    /*** SIMULATION ITERATIONS ***/

    {"iterations", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &simulation_iterations,
     0, "Number of iterations for which to run the simulation", "INT"},

    {"round_trips", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &simulation_round_trips,
     0, "Number of round trips (pessimistic samples) to run the simulation", "INT"},

    /*** MONTE CARLO OPTIMIZATION PARAMETERS ***/

    {"neighbor_scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &neighbor_scale, 0, "Ratio of neighbor sphere radius to interaction scale "
     "times ball radius. Drastically reduces collision detections","DOUBLE"},
    {"translation_scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &sw.translation_scale, 0, "Standard deviation for translations of balls, "
     "relative to ball radius", "DOUBLE"},
    {"acceptance_goal", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &acceptance_goal, 0, "Goal to set the acceptance rate", "DOUBLE"},

    /*** PARAMETERS DETERMINING OUTPUT FILE DIRECTORY AND NAMES ***/

    {"data_dir", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &data_dir, 0,
     "Directory in which to save data", "data_dir"},
    {"filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     "Base of output file names", "STRING"},
    {"filename_suffix", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
     &filename_suffix, 0, "Output file name suffix", "STRING"},

    /*** OUTPUT DATA PARAMETERS ***/

    {"de_g", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &de_g, 0,
     "Resolution of distribution functions", "DOUBLE"},
    {"de_density", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &de_density, 0, "Resolution of density file", "DOUBLE"},
    {"max_rdf_radius", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &max_rdf_radius, 0, "Set maximum radius for RDF data collection", "DOUBLE"},

    /*** HISTOGRAM METHOD OPTIONS ***/

    {"nw", '\0', POPT_ARG_NONE, &no_weights, 0, "Don't use weighing method "
     "to get better statistics on low entropy states", "BOOLEAN"},
    {"kT", '\0', POPT_ARG_DOUBLE, &fix_kT, 0, "Use a fixed temperature of kT"
     " rather than adjusted weights", "DOUBLE"},
    {"tmi", '\0', POPT_ARG_NONE, &tmi, 0,
     "Use transition matrix initialization", "BOOLEAN"},
    {"tmi-version", '\0', POPT_ARG_INT, &tmi_version, 0,
     "Use tmi version", "INT"},
    {"toe", '\0', POPT_ARG_NONE, &toe, 0,
     "Use transition optimized ensemble", "BOOLEAN"},
    {"tmmc", '\0', POPT_ARG_NONE, &tmmc, 0,
     "Use transition matrix monte carlo", "BOOLEAN"},
    {"oetmmc", '\0', POPT_ARG_NONE, &oetmmc, 0,
     "Use optimized-ensemble transition matrix monte carlo", "BOOLEAN"},
    {"wang_landau", '\0', POPT_ARG_NONE, &wang_landau, 0,
     "Use Wang-Landau histogram method", "BOOLEAN"},
    {"wltmmc", '\0', POPT_ARG_NONE, &wltmmc, 0,
     "Use Wang-Landau TMMC method", "BOOLEAN"},
    {"vanilla_wang_landau", '\0', POPT_ARG_NONE, &vanilla_wang_landau, 0,
     "Use Wang-Landau histogram method with vanilla settings", "BOOLEAN"},
    {"simple_flat", '\0', POPT_ARG_NONE, &simple_flat, 0,
     "Use simple flat histogram method", "BOOLEAN"},
    {"optimized_ensemble", '\0', POPT_ARG_NONE, &optimized_ensemble, 0,
     "Use a optimized ensemble weight histogram method", "BOOLEAN"},
    {"transition_override", '\0', POPT_ARG_NONE, &transition_override, 0,
     "Override initialized weights with weights generated from the transition matrix",
     "BOOLEAN"},

    {"transitions_input_filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
     &transitions_input_filename, 0, "File from which to read in transition matrix, "
     "which will be used to generate a weight array in place of initialization", "STRING"},
    {"golden", '\0', POPT_ARG_NONE, &golden, 0, "Run gold standard calculation", "BOOLEAN"},

    /*** HISTOGRAM METHOD PARAMETERS ***/

    {"wl_factor", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &wl_factor,
     0, "Initial value of Wang-Landau factor", "DOUBLE"},
    {"wl_fmod", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &wl_fmod, 0,
     "Wang-Landau factor modifiction parameter", "DOUBLE"},
    {"wl_threshold", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &wl_threshold, 0, "Threhold for normalized standard deviation in "
     "energy histogram at which to adjust Wang-Landau factor", "DOUBLE"},
    {"wl_cutoff", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &wl_cutoff, 0, "Cutoff for Wang-Landau factor", "DOUBLE"},
    {"oe_update_factor", '\0', POPT_ARG_INT, &oe_update_factor, 0,
     "Update scaling for the optimized ensemble method", "INT"},
    {"flat_update_factor", '\0', POPT_ARG_INT, &flat_update_factor, 0,
     "Update scaling for the simple flat method", "INT"},
    {"min_important_energy", '\0', POPT_ARG_INT, &sw.min_important_energy, 0,
     "Fix a minimum important energy at a given value", "INT"},
    {"max_entropy_energy", '\0', POPT_ARG_INT, &sw.max_entropy_state, 0,
     "Set the maximum-entropy energy", "INT"},

    /*** END CONDITION PARAMETERS ***/

    {"optimistic_sampling", '\0', POPT_ARG_VAL, &optimistic_sampling, 1,
     "Sample optimistically", 0},
    {"pessimistic_sampling", '\0', POPT_ARG_VAL, &optimistic_sampling, 0,
     "Sample pessimistically", 0},
    {"min_samples", '\0', POPT_ARG_INT, &sw.min_samples, 0,
     "Number of times to sample mininum energy", "INT"},
    {"sample_error", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &sw.sample_error, 0,
     "Fractional sample error to acquired to exit initialization", "DOUBLE"},
    {"flatness", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &sw.flatness, 0,
     "Maximum allowable proportional deviation from mean histogram value after "
     "initialization", "DOUBLE"},
    {"init_iters", '\0', POPT_ARG_INT, &sw.init_iters, 0,
     "Number of iterations for initialization", "INT"},

    {"min_T", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &sw.min_T, 0, "The minimum temperature that we care about", "DOUBLE"},
    {"dos_precision", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &sw.fractional_dos_precision, 0,
     "Precision factor for computing weights from the transition matrix", "DOUBLE"},

    /*** TESTING AND DEBUGGING OPTIONS ***/

    {"R", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &R, 0, "Ball radius (for testing purposes; should always be 1)", "DOUBLE"},
    {"seed", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0,
     "Seed for the random number generator", "INT"},
    {"debug", '\0', POPT_ARG_NONE, &debug, 0, "Debug mode", "BOOLEAN"},

    POPT_AUTOHELP
    POPT_TABLEEND
  };
  optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTION...]\nRequired arguments: number of balls (N), "
                         "and either filling fraction (ff) "
                         "or cell dimensions (lenx, leny, lenz).");

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
  printf("\n");
  if (debug) printf("DEBUG MODE IS ENABLED!\n");
  printf("------------------------------------------------------------------\n\n");

  // ----------------------------------------------------------------------------
  // Verify we have valid and reasonable arguments, and set secondary parameters
  // ----------------------------------------------------------------------------

  // Make sure we have a valid square-well fluid
  if(sw.well_width < 1){
    printf("Well width should be >= 1 (but is %g).\n",sw.well_width);
    return 254;
  }

  // Check code compatability with the requested number of dimensions with walls
  if(sw.walls >= 2){
    printf("Code cannot currently handle walls in more than one dimension.\n");
    return 254;
  }
  if(sw.walls > 3){
    printf("You cannot have walls in more than three dimensions.\n");
    return 254;
  }

  const bool reading_in_transition_matrix = (strcmp(transitions_input_filename,"none") != 0);

  // Check that only one histogram method is used
  if(bool(no_weights) + bool(simple_flat) + bool(wang_landau)
     + vanilla_wang_landau + tmi + toe + wltmmc + tmmc + oetmmc + (fix_kT != 0)
     + reading_in_transition_matrix + golden != 1){
    printf("Exactly one histogram method must be selected!\n");
    return 254;
  }
  // Check that no more than one secondary method is used
  if(optimized_ensemble + transition_override > 1){
    printf("Cannot use more than one secondary histogram method!\n");
    return 254;
  }

  // use tmmc for golden calculations
  if (golden) tmi = true;

  // Check that we are only using one end condition
  if(sw.min_samples && sw.sample_error && sw.flatness){
    printf("Can only use one end condition!\n");
    return 157;
  }

  /* If we are going to optimized the ensemble after initializing via some other method,
     initialize "half way" each time */
  if(optimized_ensemble && sw.flatness){
    printf("It does not make sense to optimize the ensemble with a "
           "flat histogram end condition!\n");
    return 175;
  }

  if (sw.len[x]<0 || sw.len[y]<0 || sw.len[x]<0){
    printf("Cell dimensions cannot be negative.\n");
    return 124;
  }

  if (sw.len[x] + sw.len[y] + sw.len[z] == 0) {
    const double volume = 4*M_PI/3*R*R*R*sw.N/sw.filling_fraction;
    const double min_cell_width = 2*sqrt(2)*R; // minimum cell width
    const int numcells = (sw.N+3)/4; // number of unit cells we need
    const int max_cubic_width
      = pow(volume/min_cell_width/min_cell_width/min_cell_width, 1.0/3);
    if (max_cubic_width*max_cubic_width*max_cubic_width >= numcells) {
      // We can get away with a cubic cell, so let's do so.
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
  } else if (sw.len[x]==0 || sw.len[y]==0 || sw.len[x]==0){
    printf("Must specify either filling fraction or ALL cell dimensions.\n");
    return 44;
  } else {
    sw.filling_fraction = (4*M_PI/3*R*R*R*sw.N)/(sw.len[x]*sw.len[y]*sw.len[z]);
  }

  printf("\nSetting cell dimensions to (%g, %g, %g).\n",
         sw.len[x], sw.len[y], sw.len[z]);

  // Compute our actual filling fraction, eta
  const double eta = (double)sw.N*4.0/3.0*M_PI*R*R*R/(sw.len[x]*sw.len[y]*sw.len[z]);
  if (eta > 1) {
    fprintf(stderr, "\nYou're trying to cram too many balls into the cell. "
            "They will never fit. Filling fraction: %g\n", eta);
    return 7;
  }

  if (sw.translation_scale == 0) {
    const double dist_in_well = 2*R*(sw.well_width-1);
    if (dist_in_well > 0.2*R) {
      sw.translation_scale = 0.1*R;
    } else {
      sw.translation_scale = dist_in_well/2;
    }
  }

  if (sw.N <= 0 || simulation_iterations < 0 ||
      simulation_round_trips < 0 || R <= 0 ||
      neighbor_scale <= 0 || sw.translation_scale < 0 ||
      sw.len[x] < 0 || sw.len[y] < 0 || sw.len[z] < 0) {
    fprintf(stderr, "\nAll parameters must be positive.\n");
    return 1;
  }

  // Choose necessary but unspecified parameters
  if(tmmc || oetmmc | tmi | toe){
    sw.sim_dos_type = transition_dos;
  } else {
    sw.sim_dos_type = histogram_dos;
  }

  /* set end condition defaults */
  if(sw.end_condition == none && !fix_kT && !no_weights && !reading_in_transition_matrix){
    // This is the default default, which may be overridden below by a
    // different default for given algorithms.
    sw.end_condition = init_iter_limit;
  }

  // set unspecified end condition parameters
  if (sw.end_condition == optimistic_min_samples && !sw.min_samples) {
    sw.min_samples = default_optimistic_min_samples;
    // different default for tmmc methods because they keep
    // accumulating statistics without doing resets.
    if (tmmc || oetmmc) sw.min_samples = tmmc_min_samples;
    printf("Defaulting min_samples to %d\n", sw.min_samples);
  } else if (sw.end_condition == pessimistic_min_samples && !sw.min_samples) {
    sw.min_samples = default_pessimistic_min_samples;
  }
  else if((sw.end_condition == optimistic_sample_error ||
           sw.end_condition == pessimistic_sample_error) &&
          !sw.sample_error){
    sw.sample_error = optimistic_sampling ?
      default_optimistic_sample_error : default_pessimistic_sample_error;
  } else if(sw.end_condition == flat_histogram && !sw.flatness){
    sw.flatness = default_flatness;
  } else if(sw.end_condition == init_iter_limit && !sw.init_iters){
    sw.init_iters = simulation_iterations;
  }

  // Set end condition text
  char *end_condition_text = new char[1024];

  if(sw.min_samples && optimistic_sampling){
    sw.end_condition = optimistic_min_samples;
    sprintf(end_condition_text,"optimistic_min_samples");
  } else if(sw.min_samples && !optimistic_sampling) {
    sw.end_condition = pessimistic_min_samples;
    sprintf(end_condition_text,"pessimmistic_min_samples");
  } else if(sw.sample_error && optimistic_sampling) {
    sw.end_condition = optimistic_sample_error;
    sprintf(end_condition_text,"optimistic_sample_error");
  } else if(sw.sample_error && !optimistic_sampling) {
    sw.end_condition = pessimistic_sample_error;
    sprintf(end_condition_text,"pessimistic_sample_error");
  } else if(sw.flatness){
    sw.end_condition = flat_histogram;
    sprintf(end_condition_text,"flat_histogram");
  } else if(sw.init_iters){
    sw.end_condition = init_iter_limit;
    sprintf(end_condition_text,"init_iter_limit");
  } else {
    sw.end_condition = none;
    sprintf(end_condition_text,"none");
  }

  // Set default data directory
  if (strcmp(data_dir,"none") == 0){
    sprintf(data_dir,"%s/s%03ld",default_data_dir,seed);
    printf("\nUsing default data directory: [deft]/%s\n",data_dir);
  }

  // If a filename was not selected, make a default
  if (strcmp(filename, "none") == 0) {
    char *method_tag = new char[200];
    char *wall_tag = new char[100];
    if(sw.walls == 0) sprintf(wall_tag,"periodic");
    else if(sw.walls == 1) sprintf(wall_tag,"wall");
    else if(sw.walls == 2) sprintf(wall_tag,"tube");
    else if(sw.walls == 3) sprintf(wall_tag,"box");
    if (fix_kT) {
      sprintf(method_tag, "-kT%g", fix_kT);
    } else if (no_weights) {
      sprintf(method_tag, "-nw");
    } else if (simple_flat) {
      sprintf(method_tag, "-simple_flat");
    } else if (tmi) {
      if (tmi_version == 1) sprintf(method_tag, "-tmi");
      else sprintf(method_tag, "-tmi%d", tmi_version);
    } else if (toe) {
      if (tmi_version == 1) sprintf(method_tag, "-toe");
      else sprintf(method_tag, "-toe%d", tmi_version);
    } else if (tmmc) {
      sprintf(method_tag, "-tmmc");}
      else if (wltmmc) {
      sprintf(method_tag, "-wltmmc");
    } else if (oetmmc) {
      sprintf(method_tag, "-oetmmc");
    } else if (wang_landau) {
      sprintf(method_tag, "-wang_landau");
    } else if (vanilla_wang_landau) {
      sprintf(method_tag, "-vanilla_wang_landau");
    } else if (reading_in_transition_matrix) {
      sprintf(method_tag, "-cfw");
    } else {
      printf("We could not identify a method for a method tag.\n");
      return 104;
    }

    if(optimized_ensemble){
      sprintf(method_tag, "%s_oe", method_tag);
    } else if(transition_override){
      sprintf(method_tag, "%s_to", method_tag);
    }

    if (golden) sprintf(method_tag, "%s-golden", method_tag);

    sprintf(filename, "%s-ww%04.2f-ff%04.2f-N%i%s",
            wall_tag, sw.well_width, eta, sw.N, method_tag);
    printf("Using default file name: ");
    delete[] method_tag;
    delete[] wall_tag;
  }
  else
    printf("\nUsing given file name: ");

  // If a filename suffix was specified, add it
  if (strcmp(filename_suffix, "none") != 0)
    sprintf(filename, "%s-%s", filename, filename_suffix);
  printf("%s\n\n",filename);

  // Initialize the random number generator with our seed
  random::seed(seed);

  // ----------------------------------------------------------------------------
  // Define sw_simulation variables
  // ----------------------------------------------------------------------------

  sw.balls = new ball[sw.N];
  sw.iteration = 0;
  sw.max_entropy_state = 0;
  sw.min_energy_state = 0;

  // initialize ball radii
  for(int i = 0; i < sw.N; i++)
    sw.balls[i].R = R;

  // scale distances by ball radius
  sw.translation_scale *= R;

  // neighbor radius should scale with radius and interaction scale
  sw.neighbor_R = neighbor_scale*R*sw.well_width;

  // Find the upper limit to the maximum number of neighbors a ball
  // could have.  Then double it, because we might not yet realize
  // that some of our neighbors move out of range, and then find we
  // have even more neighbors when another ball comes into range.
  sw.max_neighbors = 2*max_balls_within(2+neighbor_scale*sw.well_width);

  // Energy histogram and weights
  sw.interaction_distance = 2*R*sw.well_width;
  sw.energy_levels = sw.N*max_balls_within(sw.interaction_distance*1.1)/2 + 1; // the 1.1 is a fudge factor
  printf("Energy levels are %d\n",sw.energy_levels);
  sw.energy_histogram = new long[sw.energy_levels]();
  sw.ln_energy_weights = new double[sw.energy_levels]();

  // Observed and sampled energies
  sw.optimistic_samples = new long[sw.energy_levels]();
  sw.pessimistic_samples = new long[sw.energy_levels]();
  sw.pessimistic_observation = new bool[sw.energy_levels]();

  // Transitions from one energy to another
  /* We are more conservative when computing the "biggest energy
     transition" than when we set energy_levels, since a single
     instance of a ball interacting with more than "biggest energy
     transition" spheres could cause a crash, while energy_levels
     would require that the average number of interactions is greater
     than this, which probably is not possible (for a significant
     number of spheres). */
  sw.biggest_energy_transition = max_balls_within(sw.interaction_distance + 1);
  sw.transitions_table =
    new long[sw.energy_levels*(2*sw.biggest_energy_transition+1)]();

  // Walker histograms
  sw.walkers_up = new long[sw.energy_levels]();

  // ----------------------------------------------------------------------------
  // Define data arrays
  // ----------------------------------------------------------------------------

  // Radial distribution function (RDF) histogram
  long *g_energy_histogram = new long[sw.energy_levels]();
  const int g_bins = round(min(min(min(sw.len[y],sw.len[z]),sw.len[x])/2,max_rdf_radius)
                           / de_g);
  long **g_histogram = new long*[sw.energy_levels];
  for(int i = 0; i < sw.energy_levels; i++)
    g_histogram[i] = new long[g_bins]();

  // Density histogram
  const int density_bins = round(sw.len[wall_dim]/de_density);
  long **density_histogram = new long*[sw.energy_levels];
  for(int i = 0; i < sw.energy_levels; i++)
    density_histogram[i] = new long[density_bins]();

  /* Set the transitions filename so we can use it during initialization! */
  sw.transitions_filename = new char[1024];
  sprintf((char *)sw.transitions_filename, "%s/%s-transitions.dat", data_dir, filename);

  sw.transitions_movie_filename_format = new char[1024];
  sprintf((char *)sw.transitions_movie_filename_format, "%s/%s-movie", data_dir, filename);
  mkdir(sw.transitions_movie_filename_format, 0777);
  sprintf((char *)sw.transitions_movie_filename_format,
          "%s/%s-movie/%%06d-transitions.dat", data_dir, filename);

  sw.dos_movie_filename_format = new char[1024];
  sprintf((char *)sw.dos_movie_filename_format,
          "%s/%s-movie/%%06d-lndos.dat", data_dir, filename);

  sw.lnw_movie_filename_format = new char[1024];
  sprintf((char *)sw.lnw_movie_filename_format,
          "%s/%s-movie/%%06d-lnw.dat", data_dir, filename);

  char *histogram_movie_filename_format = new char[1024];
  sprintf(histogram_movie_filename_format,
          "%s/%s-movie/%%06d-E.dat", data_dir, filename);

  // ----------------------------------------------------------------------------
  // Set up the initial grid of balls
  // ----------------------------------------------------------------------------

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

  // If we made our cells too small, return with error
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

  // ----------------------------------------------------------------------------
  // Make sure initial placement is valid
  // ----------------------------------------------------------------------------

  for(int i = 0; i < sw.N; i++) {
    for(int j = 0; j < i; j++) {
      if (overlap(sw.balls[i], sw.balls[j], sw.len, sw.walls)) {
        print_bad(sw.balls, sw.N, sw.len, sw.walls);
        printf("Error in initial placement: balls are overlapping.\n");
        printf("min_cell_width = %g\n", min_cell_width);
        printf("cells = %d x %d x %d\n", cells[0], cells[1], cells[2]);
        printf("len = %g x %g x %g\n", sw.len[0], sw.len[1], sw.len[2]);
        return 253;
      }
    }
  }

  fflush(stdout);

  // ----------------------------------------------------------------------------
  // Initialization of cell and histogram methods
  // ----------------------------------------------------------------------------

  sw.energy =
    count_all_interactions(sw.balls, sw.N, sw.interaction_distance, sw.len,
                           sw.walls, sw.sticky_wall);

  // First, let us figure out what the max entropy point is (and move to it)
  sw.max_entropy_state = sw.initialize_max_entropy(acceptance_goal);
  sw.reset_histograms();
  sw.iteration = 0;

  // Now let's initialize our weight array
  if (toe || tmi || golden || wang_landau) {
    sprintf(transitions_input_filename, "%s/%s-transitions.dat", data_dir, filename);

    FILE *transitions_infile = fopen(transitions_input_filename,"r");
    if (transitions_infile != NULL) {
      fclose(transitions_infile);
      sw.initialize_transitions_file(transitions_input_filename);
      seed = random::seed_randomly();
      printf("Initializing from transitions file '%s' and using random seed %ld\n",
             transitions_input_filename, seed);
    } else {
      printf("NOT initializing from transitions file '%s', since it doesn't seem to exist\n",
             transitions_input_filename);
    }
  }

  if (reading_in_transition_matrix){
    sw.initialize_transitions_file(transitions_input_filename);
  } else if (fix_kT) {
    sw.initialize_canonical(fix_kT);
  } else if (wang_landau || vanilla_wang_landau) {
    sw.wl_factor = wl_factor;
    sw.initialize_wang_landau(wl_fmod, wl_threshold, wl_cutoff,
                              vanilla_wang_landau);
  } else if (wltmmc) {
    sw.wl_factor = wl_factor;
    sw.initialize_wltmmc(wl_fmod, wl_threshold, wl_cutoff);
  } else if (simple_flat) {
    sw.initialize_simple_flat(flat_update_factor);
  } else if (tmi) {
    sw.initialize_tmi();
  } else if (toe) {
    sw.initialize_toe();
  } else if (tmmc) {
    sw.use_tmmc = true;
    sw.initialize_transitions();
  } else if (oetmmc) {
    sw.initialize_transitions();
    sw.optimize_weights_using_transitions(tmi_version);
  }

  // If we wish to optimize the ensemble or set transition matrix weights, do so
  if (optimized_ensemble) {
    printf("\nOptimizing the ensemble!\n");
    // We need to know the minimum important energy to get optimized_ensemble right.
    if(!vanilla_wang_landau) sw.set_min_important_energy(); // because vanilla fixes min_e
    delete[] sw.compute_walker_density_using_transitions(); // to print the sample rate
    sw.reset_histograms();
    sw.iteration = 0;

    // a guess for the number of iterations for which to initially run
    //   optimized ensemble initialization
    int first_update_iterations = sw.N*sw.energy_levels;

    sw.initialize_optimized_ensemble(first_update_iterations, oe_update_factor);
    delete[] sw.compute_walker_density_using_transitions(); // to print the sample rate
  } else if(transition_override){
    printf("\nOptimizing the weight array using the transition matrix!\n");
    sw.set_min_important_energy();
    sw.optimize_weights_using_transitions(tmi_version);
    took("Optimizing with tmmc");
  }
  sw.flush_weight_array();

  took("Actual initialization");

  if (!fix_kT && sw.min_T > 0) {
    /* Force canonical weights at low energies */
    sw.set_min_important_energy();
    printf("Using canonical energies below %d (E/N = %g)\n",
           sw.min_important_energy, -sw.min_important_energy/double(sw.N));
    sw.initialize_canonical(sw.min_T,sw.min_important_energy);
  }

  if (false) {
    int E1 = sw.max_entropy_state;
    int E2 = sw.min_energy_state;
    switch (sw.N) {
    case 20:
      E2 = 95;
      break;
    }
    printf("Round trip should take %g and %g moves going down and up from %d to %d.\n",
           sw.estimate_trip_time(E1, E2), sw.estimate_trip_time(E2, E1), E1, E2);
  }

  double fractional_sample_error = sw.fractional_sample_error(sw.min_T,optimistic_sampling);

  // ----------------------------------------------------------------------------
  // Generate save file info
  // ----------------------------------------------------------------------------

  mkdir(data_dir, 0777); // create save directory

  char *e_fname = new char[1024];
  sprintf(e_fname, "%s/%s-E.dat", data_dir, filename);

  char *w_fname = new char[1024];
  sprintf(w_fname, "%s/%s-lnw.dat", data_dir, filename);

  char *os_fname = new char[1024];
  sprintf(os_fname, "%s/%s-os.dat", data_dir, filename);

  char *ps_fname = new char[1024];
  sprintf(ps_fname, "%s/%s-ps.dat", data_dir, filename);

  char *density_fname = new char[1024];
  sprintf(density_fname, "%s/%s-density.dat", data_dir, filename);

  char *headerinfo = new char[4096];
  sprintf(headerinfo,
          "# version: %s\n"
          "# well_width: %g\n"
          "# ff: %.16g\n"
          "# N: %i\n"
          "# walls: %i\n"
          "# cell dimensions: (%g, %g, %g)\n"
          "# seed: %lu\n"
          "# de_g: %g\n"
          "# de_density: %g\n"
          "# translation_scale: %g\n"
          "# neighbor_scale: %g\n"
          "# energy_levels: %i\n"
          "# min_T: %g\n"
          "# fractional_sample_error after initialization: %g\n"
          "# max_entropy_state: %d\n"
          "# min_important_energy after initialization: %i\n\n",
          version_identifier(),
          sw.well_width, sw.filling_fraction, sw.N, sw.walls, sw.len[0], sw.len[1],
          sw.len[2], seed, de_g, de_density, sw.translation_scale, neighbor_scale,
          sw.energy_levels, sw.min_T, fractional_sample_error,
          sw.max_entropy_state,
          sw.set_min_important_energy());


  if(reading_in_transition_matrix){
    sprintf(headerinfo, "%s# transitions_input_filename: %s\n",
            headerinfo, transitions_input_filename);
  }

  if(no_weights){
    sprintf(headerinfo, "%s# histogram method: none\n\n", headerinfo);
  } else if(fix_kT){
    sprintf(headerinfo,
            "%s# histogram method: canonical (fixed temperature)\n"
            "# kT: %g\n",
            headerinfo, fix_kT);
  } else if(wang_landau){
    sprintf(headerinfo,
            "%s# histogram method: Wang-Landau\n"
            "# wl_factor: %g\n"
            "# wl_fmod: %g\n"
            "# wl_threshold: %g\n"
            "# wl_cutoff: %g\n",
            headerinfo, wl_factor, wl_fmod, wl_threshold, wl_cutoff);
  } else if(vanilla_wang_landau){
    sprintf(headerinfo,
            "%s# histogram method: Vanilla Wang-Landau\n"
            "# wl_factor: %g\n"
            "# wl_fmod: %g\n"
            "# wl_threshold: %g\n"
            "# wl_cutoff: %g\n",
            headerinfo, wl_factor, wl_fmod, wl_threshold, wl_cutoff);
  } else if(simple_flat){
    sprintf(headerinfo,
            "%s# histogram method: simple flat\n", headerinfo);
  } else if (tmmc){
    sprintf(headerinfo,
            "%s# histogram method: tmmc\n",
            headerinfo);
  } else if (oetmmc){
    sprintf(headerinfo, "%s# histogram method: oetmmc\n", headerinfo);
  }

  if(sw.end_condition != none){
    sprintf(headerinfo, "%s# %s:", headerinfo, end_condition_text);
    if(sw.end_condition == optimistic_min_samples ||
       sw.end_condition == pessimistic_min_samples){
      sprintf(headerinfo, "%s %i\n", headerinfo, sw.min_samples);
    } else if(sw.end_condition == optimistic_sample_error ||
              sw.end_condition == pessimistic_sample_error){
      sprintf(headerinfo, "%s %g\n", headerinfo, sw.sample_error);
    } else if(sw.end_condition == init_iter_limit){
      sprintf(headerinfo, "%s %i\n", headerinfo, sw.init_iters);
    }
  }

  // ----------------------------------------------------------------------------
  // Print initialization info
  // ----------------------------------------------------------------------------

  char *countinfo = new char[4096];
  sprintf(countinfo,
          "# iterations: %li\n"
          "# working moves: %li\n"
          "# total moves: %li\n"
          "# acceptance rate: %g\n\n",
          sw.iteration, sw.moves.working, sw.moves.total,
          double(sw.moves.working)/sw.moves.total);

  // Save weights histogram
  FILE *w_out = fopen((const char *)w_fname, "w");
  if (!w_out) {
    fprintf(stderr, "Unable to create %s! %s\n", w_fname, strerror(errno));
    exit(1);
  }
  fprintf(w_out, "%s", headerinfo);
  fprintf(w_out, "%s", countinfo);
  fprintf(w_out, "# energy\tln(weight)\n");
  for(int i = 0; i < sw.energy_levels; i++)
    fprintf(w_out, "%i  %g\n",i,sw.ln_energy_weights[i]);
  fclose(w_out);

  delete[] countinfo;

  sw.reset_histograms();
  sw.iteration = 0;

  took("Finishing initialization");

  // ----------------------------------------------------------------------------
  // MAIN PROGRAM LOOP
  // ----------------------------------------------------------------------------

  clock_t output_period = CLOCKS_PER_SEC; // start at outputting every second
  // top out at one hour interval
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*30; // save at least every 30m
  clock_t last_output = clock(); // when we last output data

  while (sw.iteration <= simulation_iterations
         || sw.pessimistic_samples[sw.min_important_energy] < simulation_round_trips) {

    for(int i = 0; i < sw.N; i++) sw.move_a_ball();

    if (sw.iteration % (sw.N*sw.N) == 0) {
      assert(sw.energy ==
             count_all_interactions(sw.balls, sw.N, sw.interaction_distance, sw.len,
                                    sw.walls, sw.sticky_wall));
    }

    // ---------------------------------------------------------------
    // Add data to density and RDF histograms
    // ---------------------------------------------------------------

    // Density histogram -- only handles walls in one dimension
    if(sw.walls == 1){
      for(int i = 0; i < sw.N; i++){
        density_histogram[sw.energy]
          [int(floor(sw.balls[i].pos[wall_dim]/de_density))] ++;
      }
    }

    // RDF
    if (!sw.walls && sw.iteration % sw.N == 0) {
      g_energy_histogram[sw.energy]++;
      for(int i = 0; i < sw.N; i++){
        for(int j = 0; j < sw.N; j++){
          if(i != j){
            const vector3d r = periodic_diff(sw.balls[i].pos, sw.balls[j].pos, sw.len,
                                             sw.walls);
            const int r_i = floor(r.norm()/de_g);
            if(r_i < g_bins) g_histogram[sw.energy][r_i]++;
          }
        }
      }
    }

    // ---------------------------------------------------------------
    // Save data to files
    // ---------------------------------------------------------------

    // clock() can be expensive, so let's only check the time every so
    // often:
    const clock_t now = (sw.iteration % 10000) ? last_output : clock();
    if ((now - last_output > output_period)
        || (sw.iteration >= simulation_iterations
            && (simulation_round_trips == 0
                || sw.pessimistic_samples[sw.min_important_energy] >= simulation_round_trips))) {
      last_output = now;
      assert(last_output);
      output_period *= 2;
      if (output_period > max_output_period) output_period = max_output_period;
      long fraction_done = 100*sw.iteration/simulation_iterations;
      if (simulation_round_trips > 0) {
        long pessimistic_fraction =
          100*sw.pessimistic_samples[sw.min_important_energy]/simulation_round_trips;
        if (pessimistic_fraction < fraction_done) fraction_done = pessimistic_fraction;
      }
      const double secs_done = double(now)/CLOCKS_PER_SEC;
      const int seconds = int(secs_done) % 60;
      const int minutes = int(secs_done / 60) % 60;
      const int hours = int(secs_done / 3600) % 24;
      const int days = int(secs_done / 86400);
      printf("Saving data after %i days, %02i:%02i:%02i, %li iterations"
             " (%ld%%) complete, current energy %g.\n",
             days, hours, minutes, seconds, sw.iteration,
             fraction_done, -sw.energy/double(sw.N));
      fflush(stdout);

      char *countinfo = new char[4096];
      sprintf(countinfo,
              "# iterations: %li\n"
              "# working moves: %li\n"
              "# total moves: %li\n"
              "# acceptance rate: %g\n\n",
              sw.iteration, sw.moves.working, sw.moves.total,
              double(sw.moves.working)/sw.moves.total);

      // Save energy histogram
      FILE *e_out = fopen((const char *)e_fname, "w");
      fprintf(e_out, "%s", headerinfo);
      fprintf(e_out, "%s", countinfo);
      fprintf(e_out, "# max_entropy_state: %d\n",sw.max_entropy_state);
      fprintf(e_out, "# min_important_energy: %i\n\n",sw.min_important_energy);
      // fprintf(e_out, "# min_important_energy: %i\n\n",sw.set_min_important_energy());
      fprintf(e_out, "# energy   counts\n");
      for(int i = 0; i < sw.energy_levels; i++){
        if(sw.energy_histogram[i] != 0)
          fprintf(e_out, "%i  %ld\n",i,sw.energy_histogram[i]);
      }
      fclose(e_out);

      if (histogram_movie_filename_format && sw.transitions_movie_filename_format) {
        char *fname = new char[4096];
        sprintf(fname, histogram_movie_filename_format, sw.transitions_movie_count - 1);
        FILE *e_out = fopen(fname, "w");
        fprintf(e_out, "%s", headerinfo);
        fprintf(e_out, "%s", countinfo);
        fprintf(e_out, "# max_entropy_state: %d\n",sw.max_entropy_state);
        fprintf(e_out, "# min_important_energy: %i\n\n",sw.min_important_energy);
        // fprintf(e_out, "# min_important_energy: %i\n\n",sw.set_min_important_energy());
        fprintf(e_out, "# energy   counts\n");
        for(int i = 0; i < sw.energy_levels; i++){
          if(sw.energy_histogram[i] != 0)
            fprintf(e_out, "%i  %ld\n",i,sw.energy_histogram[i]);
        }
        fclose(e_out);
        delete[] fname;
      }

      // Save optimistic sample counts
      FILE *os_out = fopen(os_fname, "w");
      if (!os_out) {
        fprintf(stderr, "Unable to create %s!\n", os_fname);
        exit(1);
      }
      fprintf(os_out, "%s", headerinfo);
      fprintf(os_out, "%s", countinfo);
      fprintf(os_out, "# energy\tsamples\n");
      for(int i = 0; i < sw.energy_levels; i++) {
        if (sw.energy_histogram[i] != 0)
          fprintf(os_out, "%i  %li\n", i, sw.optimistic_samples[i]);
      }
      fclose(os_out);

      // Save pessimistic sample counts
      FILE *ps_out = fopen(ps_fname, "w");
      if (!ps_out) {
        fprintf(stderr, "Unable to create %s!\n", ps_fname);
        exit(1);
      }
      fprintf(ps_out, "%s", headerinfo);
      fprintf(ps_out, "%s", countinfo);
      fprintf(ps_out, "# energy\tsamples\n");
      for(int i = sw.max_entropy_state; i < sw.energy_levels; i++) {
        if (sw.energy_histogram[i] != 0)
          fprintf(ps_out, "%i  %li\n", i, sw.pessimistic_samples[i]);
      }
      fclose(ps_out);

      // Save RDF
      if(!sw.walls){
        char *g_fname = new char[1024];
        sprintf(g_fname, "%s/%s-g.dat", data_dir, filename);
        FILE *g_out = fopen((const char *)g_fname, "w");
        if (!g_out) {
          printf("Unable to create file %s\n", g_fname);
          exit(1);
        }
        delete[] g_fname;

        fprintf(g_out, "%s", headerinfo);
        fprintf(g_out, "%s", countinfo);
        fprintf(g_out, "# E total_counts lnw r=%g r=%g etc\n", de_g*0.5, de_g*1.5);
        fprintf(g_out, "0\t0\t0\t");
        for(int r_i = 0; r_i < g_bins; r_i++) {
          fprintf(g_out, "%g ", de_g*(r_i+0.5));
        }
        fprintf(g_out, "\n");
        for(int i = 0; i < sw.energy_levels; i++){
          if (g_energy_histogram[i] > 0){ // if we have RDF data at this energy
            fprintf(g_out, "%d\t%ld\t%g\t",
                    -i, g_energy_histogram[i], sw.ln_energy_weights[i]);
            for(int r_i = 0; r_i < g_bins; r_i++) {
              fprintf(g_out, "%ld ", g_histogram[i][r_i]);
            }
            fprintf(g_out, "\n");
          }
        }
        fclose(g_out);
      }

      // Saving density data
      if(sw.walls){
        FILE *densityout = fopen((const char *)density_fname, "w");
        fprintf(densityout, "%s", headerinfo);
        fprintf(densityout, "%s", countinfo);
        fprintf(densityout, "\n# data table containing densities in slabs "
                "(bins) of thickness de_density away from a wall");
        fprintf(densityout, "\n# row number corresponds to energy level");

        fprintf(densityout, "# E total_counts lnw z=%g z=%g etc\n",
                de_density*0.5, de_density*1.5);
        fprintf(densityout, "0\t0\t");
        for(int x_i = 0; x_i < density_bins; x_i++) {
          fprintf(densityout, "%g ", de_density*(x_i+0.5));
        }
        fprintf(densityout, "\n");
        for(int i = 0; i < sw.energy_levels; i++){
          if (sw.energy_histogram[i]) {
            fprintf(densityout, "%d\t%g",
                    -i, sw.ln_energy_weights[i]);
            for(int x_i = 0; x_i < density_bins; x_i++) {
              fprintf(densityout, "\t%ld", density_histogram[i][x_i]);
            }
            fprintf(densityout, "\n");
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

  for (int i=0; i<sw.N; i++) {
    delete[] sw.balls[i].neighbors;
  }
  delete[] sw.balls;
  delete[] sw.ln_energy_weights;
  delete[] sw.energy_histogram;

  delete[] sw.transitions_table;

  delete[] sw.walkers_up;

  delete[] sw.optimistic_samples;
  delete[] sw.pessimistic_samples;
  delete[] sw.pessimistic_observation;

  for (int i = 0; i < sw.energy_levels; i++) {
    delete[] density_histogram[i];
    delete[] g_histogram[i];
  }
  delete[] g_histogram;
  delete[] density_histogram;
  delete[] g_energy_histogram;

  delete[] headerinfo;
  delete[] e_fname;
  delete[] w_fname;
  delete[] sw.transitions_filename;
  delete[] os_fname;
  delete[] ps_fname;
  delete[] density_fname;

  delete[] data_dir;
  delete[] filename;
  delete[] filename_suffix;

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
    bool overlaps = false;
    for (int j = 0; j < i; j++) {
      if (overlap(p[i], p[j], len, walls)) {
        overlaps = true;
        break;
      }
    }
    if (overlaps) {
      char *pos = new char[1024];
      p[i].pos.tostr(pos);
      printf("%4i: %s R: %4.2f\n", i, pos, p[i].R);
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
