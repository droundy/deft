#include <stdio.h>
#include <time.h>
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <popt.h>
#include "handymath.h"
#include "vector3d.h"
#include "Monte-Carlo/square-well.h"

#define X 0
#define Y 1
#define Z 2

// -----------------------------------------------------------------------------
// Notes on conventions and definitions used
// -----------------------------------------------------------------------------
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
// Def: If two objects, a and b, are closer than a.R + b.R + neighborR + dn,
// then they are neighbors.
//
// Neighbors are used to drastically reduce the number of collision tests needed.
//
// Def: The neighborsphere of an object, a, is the sphere within which
// everything is a neighbor of a.
// Note that this sphere has a well defined center, but it does not have
// a well defined radius unless all obects are circumscribed by spheres of
// the same radius, but this does not affect anything.


// -----------------------------------------------------------------------------
// Global Constants
// -----------------------------------------------------------------------------

// prints out a lot of extra information and performs extra computation
// should be false except when something is wrong
// NOTE: This can slow things down VERY much, depending on how much debug
// code is active
const bool debug = false;



// for debugging purposes
bool testcase = false;

// -----------------------------------------------------------------------------
// Functions
// -----------------------------------------------------------------------------

// States how long it's been since last took call.
static void took(const char *name);

// Saves the locations of all balls to a file.
inline void save_locations(const ball *p, int N, const char *fname,
                           const double len[3], const char *comment="");

// The following functions only do anything if debug is true:

// Prints the location and radius of every ball
// As well as whether any overlap or are outside the cell.
inline void print_all(const ball *p, int N, double periodic[3]);

// Same as print_all, but only prints information for one ball,
// and does not do overlap tests.
inline void print_one(const ball &a, int id, const ball *p, int N,
                      double periodic[3]);

// Only print those balls that overlap or are outside the cell
// Also prints those they overlap with
inline void print_bad(const ball *p, int N, double periodic[3]);

// Checks to make sure that every ball is his neighbor's neighbor.
inline void check_neighbor_symmetry(const ball *p, int N);

int main(int argc, const char *argv[]) {
  took("Starting program");
  // -----------------------------------------------------------------------------
  // Define "Constants" -- set from arguments then unchanged
  // -----------------------------------------------------------------------------
  double periodic[3] = {0, 0, 0};
  double walls[3] = {0, 0, 0};
  bool fake_walls = false;
  unsigned long int seed = 0;

  char *dir = new char[1024];
  sprintf(dir, "papers/square-well/figs/mc");
  char *filename = new char[1024];
  sprintf(filename, "[walls/periodic]-FF");
  char *structure = new char[1024];
  sprintf(structure, "ice");
  int N = 0;
  long iterations = 100000000000;
  long initialize_iterations = 500000;
  double acceptance_goal = .4;
  double R = 1;
  double ff = 0;
  double ratio = 1;
  double neighborR = 0.5;
  double dr = 0.01;
  double de_density = 0.01;
  double de_g = 0.01;
  double dr_g = 0;
  int totime = 0;
  bool talk = false;
  // scale and theta_scale aren't quite "constants" -- they are adjusted
  // during the initialization so that we have a reasonable acceptance rate
  double scale = 0.05;
  double theta_scale = 0.05;

  poptContext optCon;
  // ---------------------------------------------------------------------------
  // Set values from parameters
  // ---------------------------------------------------------------------------
  poptOption optionsTable[] = {
    {"N", 'N', POPT_ARG_INT, &N, 0, "Number of balls to simulate", "N"},
    {"periodx", '\0', POPT_ARG_DOUBLE, &periodic[0], 0, "Periodic in x", "lenx"},
    {"periody", '\0', POPT_ARG_DOUBLE, &periodic[1], 0, "Periodic in y", "leny"},
    {"periodz", '\0', POPT_ARG_DOUBLE, &periodic[2], 0, "Periodic in z", "lenz"},
    {"wallx", '\0', POPT_ARG_DOUBLE, &walls[0], 0, "Walls in x", "lenx"},
    {"wally", '\0', POPT_ARG_DOUBLE, &walls[1], 0, "Walls in y", "leny"},
    {"wallz", '\0', POPT_ARG_DOUBLE, &walls[2], 0, "Walls in z\n", "lenz"},
    {"fake_walls", '\0', POPT_ARG_NONE, &fake_walls, 0,
     "Will cause collisions to occur with walls based on centers only", 0},
    {"iterations", 'i', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &iterations, 0,
     "Number of iterations to run for", "iter"},
    {"initialize_iterations", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT,
     &initialize_iterations, 0,
     "Number of iterations to run the initialization for", "iter"},
    {"filename", 'f', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     "Base of filename. Many files will be generated, and will automatically contain \
information on the number and type of balls, as well as file extensions", "name"},
    {"dir", 'd', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &dir, 0,
     "Directory to save to", "dir"},
    {"ratio", 'r', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &ratio, 0, "Ratio of sides of cuboid. For a cuboid of sides AxAxB, ratio = B/A.", "ratio"},
    {"R", 'R', POPT_ARG_DOUBLE, &R, 0,
     "Size of the sphere that circumscribes each ball. Defaults to setting edge length to 1", "R"},
    {"ff", '\0', POPT_ARG_DOUBLE, &ff, 0, "Filling fraction. If specified, the cell dimensions are adjusted accordingly, without changing the shape of the cell."},
    {"neighborR", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &neighborR, 0,
     "Neighbor radius, used to drastically reduce collision detections", "neighborR"},
    {"dr", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &dr, 0,
     "Differential radius change used in pressure calculation", "dr"},
    {"de_density", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &de_density, 0,
     "Resolution of density file", "de"},
    {"de_g", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &de_g, 0,
     "Resolution of distribution functions", "de"},
    {"dr_g", '\0', POPT_ARG_DOUBLE, &dr_g, 0,
     "Radius of cylinder used in distribtution functions. Defaults to the radius of a circle inscribed on a side", "dr"},
    {"scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &scale, 0,
     "Standard deviation for translations of balls", "scale"},
    {"seed", 's', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0,
     "Seed for the random number generator", "seed"},
    {"time", 't', POPT_ARG_INT, &totime, 0,
     "Timing information will be displayed", "interval"},
    {"acceptance_goal", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &acceptance_goal, 0, "Goal to set the acceptance rate", "goal"},
    POPT_AUTOHELP
    POPT_TABLEEND
  };
  optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTION...]\nNumber of balls and \
periodicity/dimensions are not set by default and are required.\n");

  int c = 0;
  // go through arguments, set them based on optionsTable
  while((c = poptGetNextOpt(optCon)) >= 0);
  if (c < -1) {
    fprintf(stderr, "\n%s: %s\n", poptBadOption(optCon, 0), poptStrerror(c));
    return 1;
  }
  poptFreeContext(optCon);

  // ---------------------------------------------------------------------------
  // Verify we have reasonable arguments and set secondary parameters
  // ---------------------------------------------------------------------------
  const bool real_walls = !fake_walls;
  double len[3] = {periodic[0] + walls[0], periodic[1] + walls[1],
                   periodic[2] + walls[2]};

  if(ff != 0) {
    // R^3*V*N/(Lx*Ly*Lz) = ff
    // Lx*Ly*Lz = R^3*V*N/ff
    // fac^3 = R^3*V*N/ff/Lx/Ly/Lz
    const double fac = R*pow(4.0/3.0*M_PI*uipow(R,3)*N/ff/len[X]/len[Y]/len[Z], 1.0/3.0);
    for(int i = 0; i < 3; i++) {
      periodic[i] *= fac;
      walls[i] *= fac;
      len[i] *= fac;
    }
    printf("\nFilling fraction was specified,");
    printf("so setting cell dimensions to (%g, %g, %g).\n", len[X], len[Y], len[Z]);

  }
  if (N <= 0 || iterations < 0 || R <= 0 || neighborR <= 0 || dr <= 0 || scale < 0) {
    fprintf(stderr, "\nAll parameters must be positive.\n");
    return 1;
  }
  for(int i = 0; i < 3; i++) {
    const char dim = i == 0 ? 'x' : i == 1 ? 'y' : 'z';
    if (periodic[i] != 0 && walls[i] != 0) {
      fprintf(stderr, "\nThe cell can't both be periodic and have walls in the %c dimension!\n", dim);
      return 1;
    }
    if (periodic[i] + walls[i] <= 0) {
      fprintf(stderr, "\nThe cell needs to have some size in the %c dimension!\n", dim);
      return 1;
    }
  }
  const double eta = (double)N*(4.0/3.0*M_PI*uipow(R,3))*R*R*R/len[X]/len[Y]/len[Z];
  if (eta > 1) {
    fprintf(stderr, "\nYou're trying to cram too many balls into the cell. They will never fit. Filling fraction: %g\n", eta);
    return 7;
  }

  if (dr_g == 0) {
    //no value was selected
    dr_g = R/2.0;
  }
  // If a filename was not selected, make a default
  if (strcmp(filename, "[walls/periodic]-FF") == 0) {
    if (walls[0] + walls[1] + walls[2] > 0) // has walls in some dimension
      sprintf(filename, "walls-%04.2f", eta);
    else
      sprintf(filename, "periodic-%04.2f", eta);
    printf("\nNo filename selected, so using the default: %s\n", filename);
  }
  printf("----------------------------------------------------------------------\n");
  printf("Running %s with parameters:\n", argv[0]);
  int cond = 0;
  for(int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') printf("\n");
    printf("%s ", argv[i]);
  }
  printf("\n");
  if (totime > 0) printf("Timing information will be displayed.\n");
  if (debug) printf("DEBUG MODE IS ENABLED!\n");
  else printf("Debug mode disabled\n");
  printf("----------------------------------------------------------------------\n\n");

  // ---------------------------------------------------------------------------
  // Define variables
  // ---------------------------------------------------------------------------
  const int density_bins = round((len[X] + len[Y] + len[Z])/de_density);
  long *density_histogram = new long[density_bins]();

  const double min_len = min(len[X], min(len[Y], len[Z]));

  const int g_bins = round(min_len/2/de_g);
  long *g_histogram = new long[g_bins]();

  ball *balls = new ball[N];
  long dZ = 0;

  // Initialize the random number generator with our seed
  random::seed(seed);

  // ---------------------------------------------------------------------------
  // Set up the initial grid of balls
  // ---------------------------------------------------------------------------
  // Approximate the maximum number of neighbors a ball could have
  // hokey guess, but it should always be high enough
  const double neighbor_sphere_vol = 4.0/3.0*M_PI*(uipow(2.5*R+neighborR,3)-uipow(R,3));
  int max_neighbors = min(N, neighbor_sphere_vol / (4.0/3.0*M_PI*uipow(R,3)));

  for(int i = 0; i < N; i++)
    balls[i].R = R;

  // Start with some lower bound for grid size
  int nx = floor(pow(N*len[X]*len[X]/(3*len[Y]*len[Z]),1.0/3.0));
  int ny = ceil(len[Y]/len[X]*nx);
  int nz = ceil(len[Z]/len[X]*nx);
  // Increase size until grid large enough to fit all balls
  //   with face-centered cubic packing
  while(fcc_total(nx,ny,nz) < N) {
    nx++;
    int ny = ceil(len[Y]/len[X]*nx);
    int nz = ceil(len[Z]/len[X]*nx);
  }

  double xspace = len[0]/(double)nx;
  double yspace = len[1]/(double)ny;
  double zspace = len[2]/(double)nz;

  // Reserve some spots at random to be vacant
  int free_spots = fcc_total(nx,ny,nz)-N;
  bool *spot_vacant = new bool[fcc_total(nx,ny,nz)]();
  int p; // Position index
  for(int i = 0; i < free_spots; i++) {
    p = floor(random::ran()*fcc_total(nx,ny,nz)); // Pick a random position index
    if(spot_vacant[p] == false)
      spot_vacant[p] = true;
    else
      i--;
  }
  // Set all remaining balls in remaining spots, keeping track of and skipping vacancies
  int encountered_vacancies = 0;
  p = 0;
  for(int i = 0; i < N; i++) {
    if(spot_vacant[i+encountered_vacancies])
      encountered_vacancies++;
    p = i + encountered_vacancies;

    // Divide the positioning problem into three subsets of all possible positions:
    //  1) The faces of the grid cells on the entire periodic cell
    //  2) The grid cell corners (the lattice) inside the entire periodic cell
    //  3) The grid cell corners on boundary of the entire periodic cell

    // If p is larger than the total number of outer and inner corner spots on grid,
    //   modulate p down to within the number of face positions (subset 1)
    if(p > fcc_outer_corners(nx,ny,nz) + fcc_inner_corners(nx,ny,nz)) {
      p -= fcc_outer_corners(nx,ny,nz) + fcc_inner_corners(nx,ny,nz);
      // Simplify problem to asking in which grid cell to position the ball
      //   and on which face of the grid cell (m)
      int m = 0;
      while(p > fcc_faces(nx,ny,nz)/3) { // Modulate p to within the number of grid cells
        p -= fcc_faces(nx,ny,nz)/3;
        m++; // Keep track how many times it has been modulated,
             //   which uniquely determines in which face to place the ball
      }
      // Set ball position in the center of a cubic grid cell
      balls[i].pos[Z] = (floor(p/(nx*ny))+0.5)*zspace;
      p = p%(nx*ny);
      balls[i].pos[Y] = (floor(p/nx)+0.5)*yspace;
      p = p%nx;
      balls[i].pos[X] = (p-1+0.5)*xspace;
      // Correct position to be on face of grid cell
      if(m == 0) { // In an x-y face of grid cell
        balls[i].pos[Z] -= 0.5*zspace;
      }
      if(m == 1) { // In a y-z face of grid cell
        balls[i].pos[X] -= 0.5*xspace;
      }
      if(m == 2) { // In a z-x face of grid cell

        balls[i].pos[Y] -= 0.5*yspace;
      }
      // If p is larger than the number of outer corner positions on the grid,
      //  but smaller than the sum of the outer and inner corner positions,
      //  modulate p down to within the number of inner corners (subset 2)
    } else if(p > fcc_outer_corners(nx,ny,nz)) {
      p -= fcc_outer_corners(nx,ny,nz);
      // For convenience, use the number of lattice positions in each dimension
      //   instead of grid positions
      nx -= 1;
      ny -= 1;
      nz -= 1;
      // Set position on inner lattice
      balls[i].pos[Z] = ceil(p/(nx*ny))*zspace;
      p = p%(nx*ny);
      balls[i].pos[Y] = ceil(p/nx)*yspace;
      p = p%nx;
      balls[i].pos[X] = p*xspace;

      // If p is within the number of outer corner positions on the grid (subset 3),
      //   split all remaining positions into four further subsets
      //   along each plane (constant: y, x, z) on the boundary of the entire cell
      //   and the periodic cell's corner
    } else if(p > nx*(ny-1)+ny*(nz-1)) { // Place in y plane boundary
      p -= nx*(ny-1)+ny*(nz-1);
      balls[i].pos[X] = (p%nz)*xspace;
      balls[i].pos[Y] = 0;
      balls[i].pos[Z] = (p%(nx-1))*zspace;
    } else if(p > nx*(ny-1)) { // Place in x plane boundary
      p -= nx*(ny-1);
      balls[i].pos[X] = 0;
      balls[i].pos[Y] = (p%(nz-1))*yspace;
      balls[i].pos[Z] = (p%ny)*zspace;
    } else if(p > 0) { // Place in constant z plane boundary
      balls[i].pos[X] = (p%(ny-1))*xspace;
      balls[i].pos[Y] = (p%nx)*yspace;
      balls[i].pos[Z] = 0;
      // If p == 0, put it on the cell's corner
    } else {
      for(int j = 0; j < 3; j++)
        balls[i].pos[j] = 0;
    }
  }
  delete[] spot_vacant;
  took("Placement");

  // ---------------------------------------------------------------------------
  // Save the initial configuration for troubleshooting
  // ---------------------------------------------------------------------------
  int most_neighbors =
    initialize_neighbor_tables(balls, N, neighborR + 2*dr, max_neighbors, periodic);
  if (most_neighbors < 0) {
    fprintf(stderr, "The guess of %i max neighbors was too low. Exiting.\n",
            max_neighbors);
    return 1;
  }
  printf("Neighbor tables initialized.");
  printf("The most neighbors is %i, whereas the max allowed is %i.\n",
         most_neighbors, max_neighbors);

  print_all(balls, N, periodic);
  print_bad(balls, N, periodic);

  // ---------------------------------------------------------------------------
  // Make sure no balls are overlapping
  // ---------------------------------------------------------------------------
  for(int i = 0; i < N; i++) {
    if (!in_cell(balls[i], walls, real_walls)) {
      printf("Oops, this is embarassing.\n");
      printf("I seem to have placed some things outside our cell.\n");
      printf("You might want to look into that.\n");
      return 17;
    }
    for(int j=i+1; j<N; j++) {
      if (overlap(balls[i], balls[j], periodic)) {
        printf("ERROR in initial placement. We have overlaps!!!\n");
        printf("AHHHHHH I DON'T KNOW WHAT TO DO!@!!!!1111\n");
        return 19;
      }
    }
  }
  fflush(stdout);

  // ---------------------------------------------------------------------------
  // Initialization of cell
  // ---------------------------------------------------------------------------
  long totalmoves = 0, workingmoves = 0, old_totalmoves = 0, old_workingmoves = 0;
  long neighbor_updates = 0, neighbor_informs = 0;
  double avg_neighbors = 0;

  double dscale = .1;

  for(long iteration=1; iteration<=initialize_iterations; iteration++) {
    // ---------------------------------------------------------------
    // Move each ball once
    // ---------------------------------------------------------------
    for(int i = 0; i < N; i++) {
      totalmoves ++;
      int move_val = move_one_ball(i, balls, N, periodic, walls, real_walls,
                                   neighborR, scale, max_neighbors, dr);
      workingmoves += move_val & 1;
      neighbor_updates += (move_val & 2) > 0;
      neighbor_informs += (move_val & 4) > 0;
    }
    // ---------------------------------------------------------------
    // fine-tune scale so that the acceptance rate will reach the goal
    // ---------------------------------------------------------------
    if (iteration % 1000 == 0) {
      const double acceptance_rate =
        (double)(workingmoves-old_workingmoves)/(totalmoves-old_totalmoves);
      old_workingmoves = workingmoves;
      old_totalmoves = totalmoves;
      if (acceptance_rate < acceptance_goal) {
        scale /= (1+dscale);
        theta_scale /= (1+scale);
      } else {
        scale *= (1+dscale);
        theta_scale *= (1+scale);
      }
      // hokey heuristic for tuning dscale
      const double closeness = fabs(acceptance_rate - acceptance_goal)/acceptance_rate;
      if(closeness > 0.5) dscale *= 2;
      else if(closeness < dscale*2) dscale/=2;
    }
    // ---------------------------------------------------------------
    // Print out timing information if desired
    // ---------------------------------------------------------------
    if (totime > 0 && iteration % totime == 0) {
      char *iter = new char[1024];
      sprintf(iter, "%i iterations", totime);
      took(iter);
      delete[] iter;
      printf("Iteration %li, acceptance rate of %g, scale: %g.\n",
             iteration, (double)workingmoves/totalmoves, scale);
      printf("We've had %g updates per kilomove and %g informs per kilomove, for %g informs per update.\n", 1000.0*neighbor_updates/totalmoves, 1000.0*neighbor_informs/totalmoves, (double)neighbor_informs/neighbor_updates);
      const long checks_without_tables = totalmoves*N;
      int total_neighbors = 0;
      for(int i = 0; i < N; i++) {
        total_neighbors += balls[i].num_neighbors;
        most_neighbors = max(balls[i].num_neighbors, most_neighbors);
      }
      avg_neighbors = double(total_neighbors)/N;
      const long checks_with_tables = totalmoves*avg_neighbors + N*neighbor_updates;
      printf("We've done about %.3g%% of the distance calculations we would have done without tables.\n", 100.0*checks_with_tables/checks_without_tables);
      printf("The max number of neighbors is %i, whereas the most we've seen is %i.\n", max_neighbors, most_neighbors);
      printf("Neighbor radius is %g and avg. number of neighbors is %g.\n\n", neighborR, avg_neighbors);
      fflush(stdout);
    }
  }
  took("Initialization");

  // ---------------------------------------------------------------------------
  // Generate header info to put in save files
  // ---------------------------------------------------------------------------
  char *headerinfo = new char[4096];
  sprintf(headerinfo, "\
# period: (%5.2f, %5.2f, %5.2f), walls: (%5.2f, %5.2f, %5.2f), de_density: %g, de_g: %g\n\
# dr_g: %g, seed: %li, R: %f, scale: %g, theta_scale: %g, real_walls: %i\n\
# initialize_iterations: %li, initial structure: %s, neighborR: %g, dr: %g\n",
          periodic[0], periodic[1], periodic[2], walls[0], walls[1], walls[2],
          de_density, de_g, dr_g, seed, R, scale, theta_scale, real_walls,
          initialize_iterations, structure, neighborR, dr);


  // fixme: to use this again, make it so file stays open and is fflushed
  // so that the name only appears once
  // // ---------------------------------------------------------------
  // // Clear the pressure file
  // // ---------------------------------------------------------------
  // char *pressure_fname = new char[1024];
  // sprintf(pressure_fname, "%s/%s-pressure-%s-%i.dat", dir, filename, shape.name, N);
  // FILE *pressureout = fopen((const char *)pressure_fname, "w");
  // const double dV = N*shape.volume*(uipow(R+dr, 3) - uipow(R, 3));
  // fprintf(pressureout, "0 0 %g #dV\n", dV);
  // fprintf(pressureout, "# total moves    pressure     dZ\n");
  // delete[] pressure_fname;
  // fclose(pressureout);

  // ---------------------------------------------------------------------------
  // MAIN PROGRAM LOOP
  // ---------------------------------------------------------------------------
  clock_t output_period = CLOCKS_PER_SEC; // start at outputting every minute
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*60; // top out at one hour interval
  clock_t last_output = clock(); // when we last output data

  int frame = 0;
  totalmoves = 0, workingmoves = 0, old_totalmoves = 0, old_workingmoves = 0;
  neighbor_updates = 0, neighbor_informs = 0;

  for(long iteration=1; iteration<=iterations; iteration++) {
    // ---------------------------------------------------------------
    // Move each ball once
    // ---------------------------------------------------------------
    for(int i = 0; i < N; i++) {
      totalmoves ++;
      int move_val = move_one_ball(i, balls, N, periodic, walls, real_walls,
                                   neighborR, scale, max_neighbors, dr);
      workingmoves += move_val & 1;
    }
    // ---------------------------------------------------------------
    // Add data to historams
    // ---------------------------------------------------------------
    for(int i = 0; i < N; i++) {
      // Density histogram:
      const int x_i = floor(balls[i].pos[X]/de_density);
      const int y_i = floor(balls[i].pos[Y]/de_density);
      const int z_i = floor(balls[i].pos[Z]/de_density);
      density_histogram[x_i] ++;
      density_histogram[int(round(len[X]/de_density)) + y_i] ++;
      density_histogram[int(round((len[X] + len[Y])/de_density)) + z_i] ++;
    }

    // ---------------------------------------------------------------
    // Get pressure info - this is slow, might want to do it less often
    // fixme - disabled for now
    // ---------------------------------------------------------------
    // for(int i = 0; i < N; i++) {
    //   dZ += overlaps_with_any(balls[i], balls, periodic, true, dr);
    // }
    // dZ /= 2; // since we check each thing against each of its neighbors, we end up
    //          // double counting. Ideally, we would count half as much instead of dividing
    // ---------------------------------------------------------------
    // Save to file
    // ---------------------------------------------------------------
    const clock_t now = clock();
    if ((now > last_output + output_period) || iteration==iterations) {
      last_output = now;
      assert(last_output);
      if (output_period < max_output_period/2) output_period *= 2;
      else if (output_period < max_output_period) output_period = max_output_period;
      const double secs_done = double(now)/CLOCKS_PER_SEC;
      const int seconds = int(secs_done) % 60;
      const int minutes = int(secs_done / 60) % 60;
      const int hours = int(secs_done / 3600) % 24;
      const int days = int(secs_done / 86400);
      printf("Saving data after %i days, %02i:%02i:%02i, %li iterations complete.\n",
             days, hours, minutes, seconds, iteration);
      fflush(stdout);

      char *countinfo = new char[4096];
      sprintf(countinfo, "\
# iteration: %li, workingmoves: %li, totalmoves: %li, acceptance rate: %g\n",
              iteration, workingmoves, totalmoves, double(workingmoves)/totalmoves);


      // Saving density in each of the x, y, z dimensions
      char *density_fname = new char[1024];
      sprintf(density_fname, "%s/%s-density-%i.dat", dir, filename, N);
      FILE *densityout = fopen((const char *)density_fname, "w");
      delete[] density_fname;
      fprintf(densityout, "%s", headerinfo);
      fprintf(densityout, "%s", countinfo);
      fprintf(densityout, "\n#e       xdensity   ydensity   zdensity   histograms in x,y,z order\n");
      const int xbins = round(len[X]/de_density);
      const int ybins = round(len[Y]/de_density);
      const int zbins = round(len[Z]/de_density);
      int maxbin = max(max(xbins, ybins), zbins);
      for(int e_i = 0; e_i < maxbin; e_i ++) {
        const double e = (e_i + 0.5)*de_density;
        const double xshell_volume = len[Y]*len[Z]*de_density;
        const double yshell_volume = len[X]*len[Z]*de_density;
        const double zshell_volume = len[X]*len[Y]*de_density;
        // If we have a non-cubic cell, then this makes the density negative
        // in the dimensions where it does not exist.
        const long xhist = e_i>=xbins ? -totalmoves : density_histogram[e_i];
        const long yhist = e_i>=ybins ? -totalmoves : density_histogram[xbins + e_i];
        const long zhist = e_i>=zbins ? -totalmoves : density_histogram[xbins + ybins + e_i];
        const double xdensity = (double)xhist*N/totalmoves/xshell_volume;
        const double ydensity = (double)yhist*N/totalmoves/yshell_volume;
        const double zdensity = (double)zhist*N/totalmoves/zshell_volume;
        fprintf(densityout, "%6.3f   %8.5f   %8.5f   %8.5f   %li %li %li\n", e, xdensity, ydensity, zdensity, xhist, yhist, zhist);
      }
      fclose(densityout);

      // Save distribution funtions
      // fixme: assumes homogeneous for density
      char *g_fname = new char[1024];
      sprintf(g_fname, "%s/%s-g-%i.dat", dir, filename, N);
      FILE *g_out = fopen((const char *)g_fname, "w");
      delete[] g_fname;
      fprintf(g_out, "%s", headerinfo);
      fprintf(g_out, "%s", countinfo);
      fprintf(g_out, "\ne       ");
      const double density = N/len[X]/len[Y]/len[Z];
      const double vol0 = len[X]*len[Y]*len[Z];
      fprintf(g_out, "\n");
      for(int e_i = 0; e_i < g_bins; e_i++) {
        const double e = (e_i + 0.5)*de_g;
        fprintf(g_out, "%6.3f  ", e);
        const double vol1 = 4.0/3.0*M_PI*(uipow(e+de_g/2, 3) - uipow(e-de_g/2, 3));
        const double probability = (double)g_histogram[e_i]/totalmoves;
        const double n2 = probability/vol0/vol1;
        const double g = n2/sqr(density)*N*N;
          fprintf(g_out, "%8.5f \n", g);
      }
      fclose(g_out);

      delete[] countinfo;
    }
    // ---------------------------------------------------------------
    // Print out timing information if desired
    // ---------------------------------------------------------------
    if (totime > 0 && iteration % totime == 0) {
      char *iter = new char[1024];
      sprintf(iter, "%i iterations", totime);
      took(iter);
      delete[] iter;
      printf("Iteration %li, acceptance rate of %g, scale: %g.\n",
             iteration, (double)workingmoves/totalmoves, scale);
      printf("We've had %g updates per kilomove and %g informs per kilomove, for %g informs per update.\n", 1000.0*neighbor_updates/totalmoves, 1000.0*neighbor_informs/totalmoves, (double)neighbor_informs/neighbor_updates);
      const long checks_without_tables = totalmoves*N;
      int total_neighbors = 0;
      for(int i = 0; i < N; i++) {
        total_neighbors += balls[i].num_neighbors;
        most_neighbors = max(balls[i].num_neighbors, most_neighbors);
      }
      avg_neighbors = double(total_neighbors)/N;
      const long checks_with_tables = totalmoves*avg_neighbors + N*neighbor_updates;
      printf("We've done about %.3g%% of the distance calculations we would have done without tables.\n", 100.0*checks_with_tables/checks_without_tables);
      printf("The max number of neighbors is %i, whereas the most we've seen is %i.\n", max_neighbors, most_neighbors);
      printf("Neighbor radius is %g and avg. number of neighbors is %g.\n\n", neighborR, avg_neighbors);
      fflush(stdout);
    }
  }
  // ---------------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // ---------------------------------------------------------------------------
  print_bad(balls, N, periodic);

  //delete[] balls; fixme
  delete[] density_histogram;
  delete[] g_histogram;

  delete[] headerinfo;
  return 0;
}
// -----------------------------------------------------------------------------
// END OF MAIN
// -----------------------------------------------------------------------------

inline void print_all(const ball *p, int N, double periodic[3]) {
  if (debug) {
    for (int i = 0; i < N; i++) {
      char *pos = new char[1024];
      p[i].pos.tostr(pos);
      printf("%4i: R: %4.2f, %i neighbors: ", i, p[i].R, p[i].num_neighbors);
      for(int j = 0; j < min(10, p[i].num_neighbors); j++)
        printf("%i ", p[i].neighbors[j]);
      if (p[i].num_neighbors > 10)
        printf("...");
      printf("\n      pos:          %s\n", pos);
    }
    printf("\n");
    fflush(stdout);
  }
}

inline void print_one(const ball &a, int id, const ball *p, int N, double periodic[3]) {
  if (debug) {
    char *pos = new char[1024];
    a.pos.tostr(pos);
    printf("%4i: R: %4.2f, %i neighbors: ", id, a.R, a.num_neighbors);
    for(int j=0; j<min(10, a.num_neighbors); j++)
      printf("%i ", a.neighbors[j]);
    if (a.num_neighbors > 10)
      printf("...");
    printf("\n      pos:          %s\n", pos);
    for (int j=0; j<N; j++) {
      if (j != id && overlap(a, p[j], periodic)) {
        p[j].pos.tostr(pos);
        printf("\t  Overlaps with %i", j);
        printf(": %s\n", pos);
      }
    }
  }
  printf("\n");
  fflush(stdout);
}

inline void print_bad(const ball *p, int N, double periodic[3]) {
  if (debug) {
    for (int i = 0; i < N; i++) {
      bool incell = true; //in_cell(p[i], walls, real_walls); fixme
      bool overlaps = false;
      for (int j=0; j<N; j++) {
        if (j != i && overlap(p[i], p[j], periodic)) {
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
        for (int j=0; j<N; j++) {
          if (j != i && overlap(p[i], p[j], periodic)) {
            p[j].pos.tostr(pos);
            printf("\t  Overlaps with %i", j);
            printf(": %s\n", pos);
          }
        }
        delete[] pos;
      }
    }
  }
  fflush(stdout);
}

inline void check_neighbor_symmetry(const ball *p, int N) {
  if (debug) {
    for(int i = 0; i < N; i++) {
      for(int j=0; j<p[i].num_neighbors; j++) {
        const int k = p[i].neighbors[j];
        bool is_neighbor = false;
        for (int l=0; l<p[k].num_neighbors; l++) {
          if (p[k].neighbors[l] == i) {
            is_neighbor = true;
            break;
          }
        }
        if(!is_neighbor) {
          printf("NEIGHBOR TABLE ERROR: %i has %i as a neighbor, but %i does not reciprocate!!!\n", i, k, k);
          testcase = true;
        }
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
    printf("%s took %.0f minutes and %g seconds.\n", name, seconds/60, fmod(seconds,60));
  } else {
    printf("%s took %g seconds..\n", name, seconds);
  }
  fflush(stdout);
  last_time = t;
}

void save_locations(const ball *p, int N, const char *fname, const double len[3], const char *comment) {
  FILE *out = fopen((const char *)fname, "w");
  fprintf(out, "# %s\n", comment);
  fprintf(out, "%g %g %g\n", len[X], len[Y], len[Z]);
  for(int i = 0; i < N; i++) {
    fprintf(out, "%6.2f %6.2f %6.2f ", p[i].pos[X], p[i].pos[Y], p[i].pos[Z]);
    fprintf(out, "\n");
  }
  fclose(out);
}
