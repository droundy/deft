#include <stdio.h>
#include <time.h>
#include "vector3d.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "handymath.h"
#include <popt.h>

// -----------------------------------------------------------------------------
// Notes on conventions and definitions used
// -----------------------------------------------------------------------------
//
// All coordinates are (presently) cartesion.
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
// Def: Iff two objects, a and b, are closer than a.R + b.R + neighborR + dn,
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
// Structs defined
// -----------------------------------------------------------------------------

// We make one poly_shape for each shape we care about (cube, tetrahedron, etc.)
// We then point to it from our polyhedra to identify that they are that shape,
// and to be able to find vertices, etc.

// Note: faces is a list of normal vectors perpendicular to the actual faces.
// The sign of this vector is unimportant for the collision detection algorithm
// used, so any two parallel faces are only counted as one face.
struct poly_shape {
  int nvertices;
  int nfaces;
  vector3d *vertices;
  vector3d *faces;
  double volume;
  char *name;

  poly_shape(const char *set_name);

  ~poly_shape() {
    delete[] vertices;
    delete[] faces;
    delete[] name;
  }
};

// This contains all of the information about a given polyhedron.
// We create and use an array of these.
struct polyhedron {
  vector3d pos;
  rotation rot;
  double R;
  const poly_shape *mypoly;
  int *neighbors;
  int num_neighbors;
  vector3d neighbor_center;
};

// -----------------------------------------------------------------------------
// Global Constants
// -----------------------------------------------------------------------------

// prints out a lot of extra information and performs extra computation
// should be false except when something is wrong
// NOTE: This can slow things down VERY much, depending on how much debug
// code is active
const bool debug = false;

// Prints time information. Useful for testing optimizations,
// or seeing how much adding complexity affects runtime.
const bool totime = false;

// for debugging purposes
bool testcase = false;

// -----------------------------------------------------------------------------
// Functions
// -----------------------------------------------------------------------------

// Sets up the vertices, faces, etc. for each polyhedron. Only run once.
void define_shapes();

// Find's the neighbors of a by comparing a's position to the center of
// everyone else's neighborsphere, where n is the index of a in bs.
inline void update_neighbors(polyhedron &a, int n, const polyhedron *bs, int N, double neighborR, double periodic[3]);

// Removes n from the neighbor table of anyone neighboring oldp but not newp.
// Adds n to the neighbor table of anyone neighboring newp but not oldp.
// Keeps the neighbor tables sorted.
inline void inform_neighbors(const polyhedron &newp, const polyhedron &oldp, int n, polyhedron *p, int N);

// Attempts to move each polyhedron once. Returns the number of successful moves.
inline int move_all(polyhedron *p, int N, double periodic[3], double walls[3], bool real_walls, double neighborR, double scale, double theta_scale, int max_neighbors);

// Check whether two polyhedra overlap
inline bool overlap(const polyhedron &a, const polyhedron &b, double periodic[3]);

// Check whether a polyhedron overlaps with any in the array, except for
// Useful for post-move testing with a temporary polyhedron.
// If count is true, it will return the total number of overlaps, otherwise
// it returns 1 if there is at least one overlap, 0 if there are none.
inline int overlaps_with_any(const polyhedron &a, const polyhedron *bs, double periodic[3], bool count=false, double dr=0);

// Check if polyhedron is outside its allowed locations
inline bool in_cell(const polyhedron &p, double walls[3], bool real_walls);

// Return the vector pointing from a to b, accounting for periodic
// boundaries
inline vector3d periodic_diff(const vector3d &a, const vector3d  &b, double periodic[3]);

// Move v in a random direction by a distance determined by a gaussian distribution
inline vector3d move(const vector3d &v, double scale, const double len[3]);

// If v is outside the cell, and there are periodic boundary condition(s), it is
// moved (appropriately) into the cell
inline vector3d fix_periodic(vector3d newv, const double len[3]);

// States how long it's been since last took call.
static void took(const char *name);

// Saves the vertices of all polyhedra to a file.
void save_locations(const polyhedron *p, int N, const char *fname, long iteration = -1);




// The following functions only do anything if debug is true:

// Prints the locations, rotations, radii, and shapes of every polyhedron
// As well as if any overlap or are outside the cell.
inline void print_all(const polyhedron *p, int N, double periodic[3]);

// Same as print_all, but only prints information for one polyhedron,
// and does not do overlap tests.
inline void print_one(const polyhedron &a, int id, const polyhedron *p, int N, double periodic[3]);

// Only print those shapes that overlap or are outside the cell
// Also prints those they overlap with
inline void print_bad(const polyhedron *p, int N, double periodic[3]);

// Checks to make sure that every polyhedron is his neighbor's neighbor.
inline void check_neighbor_symmetry(const polyhedron *p, int N);



int main(int argc, const char *argv[]) {
  // ---------------------------------------------------------------------------
  // Set up the generic polyhedra we will use
  // ---------------------------------------------------------------------------
  const poly_shape cube("cube");
  const poly_shape tetrahedron("tetrahedron");
  const poly_shape truncated_tetrahedron("truncated tetrahedron");

  // -----------------------------------------------------------------------------
  // "Constants" -- set at runtime then unchanged
  // -----------------------------------------------------------------------------
  double periodic[3] = {0, 0, 0};
  double walls[3] = {0, 0, 0};
  bool real_walls = false;
  unsigned long int seed = 0;
  int vertex_period = 0;
  poptContext optCon;
  char *shape_name = new char[1024];
  sprintf(shape_name, "cube");
  const poly_shape *shape = &cube;

  char *dir = new char[1024];
  sprintf(dir, "papers/polyhedra/figs/mc");
  char *filename = new char[1024];
  sprintf(filename, "[walls/periodic]-FF");
  int N;
  long iterations = 100000000000;
  double R = sqrt(3.0)/2.0;
  double neighborR = 0.5;
  double dr = 0.01;
  double de_density = 0.01;
  // scale and theta_scale aren't quite "constants" -- they are adjusted
  // during the initialization so that we have a reasonable acceptance rate
  double scale = 0.05;
  double theta_scale = 0.05;
  // ---------------------------------------------------------------------------
  // Set values from parameters
  // ---------------------------------------------------------------------------
  poptOption optionsTable[] = {
    {"N", 'N', POPT_ARG_INT, &N, 0, "Number of polyhedra to simulate", "N"},
    {"periodx", '\0', POPT_ARG_DOUBLE, &periodic[0], 0, "Periodic in x", "lenx"},
    {"periody", '\0', POPT_ARG_DOUBLE, &periodic[1], 0, "Periodic in y", "leny"},
    {"periodz", '\0', POPT_ARG_DOUBLE, &periodic[2], 0, "Periodic in z", "lenz"},
    {"wallx", '\0', POPT_ARG_DOUBLE, &walls[0], 0, "Walls in x", "lenx"},
    {"wally", '\0', POPT_ARG_DOUBLE, &walls[1], 0, "Walls in y", "leny"},
    {"wallz", '\0', POPT_ARG_DOUBLE, &walls[2], 0, "Walls in z\n", "lenz"},
    {"real_walls", 'r', POPT_ARG_NONE, &real_walls, 0,
     "Will cause collisions to occur with walls based on vertex locations instead of centers", 0},
    {"iterations", 'i', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &iterations, 0,
     "Number of iterations to run for", "iter"},
    {"filename", 'f', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
     "Base of filename. Many files will be generated, and will automatically contain \
information on the number and type of polyhedra, as well as file extensions", "name"},
    {"dir", 'd', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &dir, 0,
     "Directory to save to", "dir"},
    {"shape", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &shape_name, 0,
     "Type of polyhedra to use. Can be one of [cube, tetrahedron, truncated_tetrahedron]", "shape"},
    {"R", 'R', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &R, 0,
     "Size of the sphere that circumscribes each polyhedron", "R"},
    {"neighborR", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &neighborR, 0,
     "Neighbor radius, used to drastically reduce collision detections", "neighborR"},
    {"dr", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &dr, 0,
     "Differential radius change used in pressure calculation", "dr"},
    {"de_density", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &de_density, 0,
     "Resolution of density file", "de"},
    {"scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &scale, 0,
     "Standard deviation for translations of polyhedra", "scale"},
    {"theta_scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &theta_scale, 0,
     "Standard deviation for rotations of polyhedra", "theta"},
    {"seed", 's', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0,
     "Seed for the random number generator", "seed"},
    {"save_vertices", 'v', POPT_ARG_INT, &vertex_period, 0,
     "Periodically saves the vertex locations of the polyhedra, where period is a number \
of iterations", "period"},
    POPT_AUTOHELP
    POPT_TABLEEND
  };
  optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTION...]\nNumber of polyhedra and periodicity/dimensions \
are not set by default and are required.\n");

  int c = 0;
  // go through options, sets them based on optionsTable
  while((c = poptGetNextOpt(optCon)) >= 0);
  if (c < -1) {
    fprintf(stderr, "\n%s: %s\n", poptBadOption(optCon, 0), poptStrerror(c));
    return 1;
  }
  poptFreeContext(optCon);

  // ---------------------------------------------------------------------------
  // Check to make sure we have reasonable arguments and set secondary parameters
  // ---------------------------------------------------------------------------
  if (strcmp(shape_name, "cube") == 0) shape = &cube;
  else if (strcmp(shape_name, "tetrahedron") == 0) shape = &tetrahedron;
  else if (strcmp(shape_name, "truncated_tetrahedron") == 0) shape = &truncated_tetrahedron;
  else {
    fprintf(stderr, "\nInvalid shape.\n");
    return 1;
  }
  if (N <= 0 || iterations <= 0 || R <= 0 || neighborR <= 0 || dr <= 0 || scale < 0 ||
      theta_scale < 0 || vertex_period < 0) {
    fprintf(stderr, "\nAll parameters must be positive.\n");
    return 1;
  }
  for(int i=0; i<3; i++) {
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
  const double len[3] = {periodic[0] + walls[0], periodic[1] + walls[1],
                         periodic[2] + walls[2]};
  const double eta = (double)N*shape->volume/len[0]/len[1]/len[2];
  if (eta > 1.0) {
    fprintf(stderr, "\nYou're trying to cram too many polyhedra into the cell. They will never fit.\n");
    return 1;
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
  for(int i=1; i<argc; i++) {
    printf("%s ", argv[i]);
    if(strcmp(argv[i], "--real_walls") == 0)
      cond = 1 - cond;
    if(i%2 == cond) printf("\n");
  }
  printf("\n");
  if (totime) printf("Timing information will be displayed.\n");
  if (debug) printf("DEBUG MODE IS ENABLED!\n");
  else printf("Debug mode disabled\n");
  printf("----------------------------------------------------------------------\n\n");

  // ---------------------------------------------------------------------------
  // Define variables
  // ---------------------------------------------------------------------------
  const int density_bins = round((len[0] + len[1] + len[2])/de_density);
  long *density_histogram = new long[3*density_bins]();
  polyhedron *polyhedra = new polyhedron[N];
  long neighbor_updates = 0, neighbor_informs = 0;
  long totalmoves = 0, workingmoves = 0;

  // Initialize the random number generator with our seed so we get repreducable
  // results.

  // fixme

  // ---------------------------------------------------------------------------
  // Initialize the cell
  // ---------------------------------------------------------------------------
  // Start by arranging the shapes in a grid.
  // Right now is set up to simply make them evenly spaced with no rotation.
  // It will fail with an error if any of them overlap.
  // NOTE: This only really works if the cell is at least nearly cubic,
  // and if there's enough room to fit the shapes without rotation.
  // A better version, that is shape-dependent, will be required for higher
  // filling fractions.

  // NOTE: All setup assumes that we do not have a mixture,
  // but only one polyhedron.

  // Approximate the maximum number of neighbors a polyhedron could have
  const double neighbor_sphere_vol = 4.0/3.0*M_PI*uipow(2.5*R+neighborR, 3);
  int max_neighbors = min(N, neighbor_sphere_vol / shape->volume);

  for(int i=0; i<N; i++) {
    polyhedra[i].mypoly = shape;
    polyhedra[i].R = R;
    polyhedra[i].neighbors = new int[max_neighbors]();
  }
  const int nx = ceil(pow(N*len[0]*len[0]/len[1]/len[2], 1.0/3.0));
  const int ny = nx;
  const int nz = ny;
  const double xspace = len[0]/double(nx);
  const double yspace = len[1]/double(ny);
  const double zspace = len[2]/double(nz);
  printf("Setting up grid with polyhedra: (%i, %i, %i) and space: (%g, %g, %g)\n",
                   nx, ny, nz, xspace, yspace, zspace);
  if (shape == &cube) {
    bool *spot_used = new bool[nx*ny*nz]();
    for(int i=0; i<N; i++) {
      int x, y, z;
      do {
        x = floor(random::ran()*nx);
        y = floor(random::ran()*ny);
        z = floor(random::ran()*nz);
      } while(spot_used[x*ny*nz + y*nz + z]);
      spot_used[x*ny*nz + y*nz + z] = true;

      polyhedra[i].pos[0] = (x + 0.5)*xspace;
      polyhedra[i].pos[1] = (y + 0.5)*yspace;
      polyhedra[i].pos[2] = (z + 0.5)*zspace;

      polyhedra[i].neighbor_center = polyhedra[i].pos;
    }
    delete[] spot_used;
  }
  if (shape == &tetrahedron || shape == &truncated_tetrahedron) {
    int x = 0, y = 0, z = 0, phi=0;
    bool alt = false;
    const vector3d axis(0,0,1);
    for(int i=0; i<N; i++) {
      polyhedra[i].pos[0] = (x + 0.5)*xspace;
      polyhedra[i].pos[1] = (y + 0.5)*yspace + alt*R*(1.0/3.0);
      polyhedra[i].pos[2] = (z + 0.5)*zspace;
      polyhedra[i].rot = rotation(phi, axis);

      polyhedra[i].neighbor_center = polyhedra[i].pos;

      x ++;
      alt = !alt;
      phi = M_PI*alt;

      if (x >= nx) {
        x = 0; y ++;
        if (y >= ny) {
          y = 0; z ++;
        }
      }
    }
  }

  // Now let's do the initial population of the neighbor tables
  for(int i=0; i<N; i++) {
    polyhedra[i].num_neighbors = 0;
    for(int j=0; j<N; j++) {
      const bool is_neighbor = (i != j) &&
        (periodic_diff(polyhedra[i].pos, polyhedra[j].pos, periodic).normsquared() <
         uipow(polyhedra[i].R + polyhedra[j].R + neighborR, 2));
      if (is_neighbor) {
        const int index = polyhedra[i].num_neighbors;
        polyhedra[i].num_neighbors ++;
        polyhedra[i].neighbors[index] = j;
      }
    }
    if (debug && false)
      printf("Polyhedron %i has %i neighbors.\n", i, polyhedra[i].num_neighbors);
  }
  // print_all(polyhedra, N, periodic);
  print_bad(polyhedra, N, periodic);
  printf("\n\n");
  // Make sure none are overlapping:
  for(int i=0; i<N; i++) {
    if (!in_cell(polyhedra[i], walls, real_walls)) {
      printf("Oops, this is embarassing. I seem to have placed some things outside our cell.\n");
      printf("You might want to look into that.\n");
      char *vertices_fname = new char[1024];
      sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape->name, N, -1);
      save_locations(polyhedra, N, vertices_fname);
      delete[] vertices_fname;
      return 17;
    }
    for(int j=i+1; j<N; j++) {
      if (overlap(polyhedra[i], polyhedra[j], periodic)) {
        printf("ERROR in intial placement. We have overlaps!!!\n");
        printf("AHHHHHH I DON'T KNOW WHAT TO DO!@!!!!1111\n");
        char *vertices_fname = new char[1024];
        sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape->name, N, -1);
        save_locations(polyhedra, N, vertices_fname);
        delete[] vertices_fname;
        return 19;
      }
    }
  }
  fflush(stdout);

  double avg_neighbors = 0;
  int most_neighbors = 0;

  // ---------------------------------------------------------------------------
  // Initialization of cell
  // ---------------------------------------------------------------------------
  printf("---------------------------------------------------------------\n");
  printf("Beginning the initialization.\n");
  printf("---------------------------------------------------------------\n");
  long iteration = 1;
  long oldtotalmoves = 0, oldworkingmoves = 0;
  bool initialization = true;
  while (initialization) {
    if (debug) {
      print_bad(polyhedra, N, periodic);
    }
    // ---------------------------------------------------------------
    // Move each polyhedron once
    // ---------------------------------------------------------------
    totalmoves += N;
    workingmoves += move_all(polyhedra, N, periodic, walls, real_walls,
                             neighborR, scale, theta_scale, max_neighbors);

    if (debug && testcase) {
      printf("totalmoves = %li.\n", totalmoves);
      return 22;
    }
    // ---------------------------------------------------------------
    // Print timing information if desired
    // ---------------------------------------------------------------
    if (totime) {
      const int numtotime = 1000;
      if (iteration % numtotime == 0) {
        char *iter = new char[1024];
        sprintf(iter, "%i iterations (current iteration: %li)", numtotime, iteration);
        took(iter);
        printf("We've had %g updates per kilomove and %g informs per kilomove, for %g informs per update.\n", 1000.0*neighbor_updates/totalmoves, 1000.0*neighbor_informs/totalmoves, (double)neighbor_informs/neighbor_updates);
        const long checks_without_tables = totalmoves*N;
        int total_neighbors = 0;
        for(int i=0; i<N; i++) {
          total_neighbors += polyhedra[i].num_neighbors;
          most_neighbors = max(polyhedra[i].num_neighbors, most_neighbors);
        }
        avg_neighbors = double(total_neighbors)/N;
        const long checks_with_tables = totalmoves*avg_neighbors + N*neighbor_updates;
        printf("We've done about %.3g%% of the distance calculations we would have done without tables.\n", 100.0*checks_with_tables/checks_without_tables);
        printf("The max number of neighbors is %i, whereas the most we've seen is %i.\n", max_neighbors, most_neighbors);
        printf("Neighbor radius is %g and avg. number of neighbors is %g.\n\n", neighborR, avg_neighbors);
        delete[] iter;
        fflush(stdout);
      }
    }
    // ---------------------------------------------------------------
    // Figure out if we're done with initialization.
    // Since more start at lower z-values, once there are more at
    // higher values, we should be readyish.
    // ---------------------------------------------------------------
    if (iteration % 1000 == 0) {
      int count1 = 0, count2 = 0;
      for(int i=0; i<N; i++) {
        if (polyhedra[i].pos[2] < len[2]/2.0)
          count1 ++;
        else
          count2 ++;
      }
      if (count1 > count2) {
          printf("Not ready to start storing data yet - c1: %i, c2: %i, iteration: %li, acceptance rate: %4.2f\n", count1, count2, iteration, (double)workingmoves/totalmoves);
      } else {
        printf("\nTime to start data collection! c1: %i, c2: %i, iteration: %li\n", count1, count2, iteration);
        took("Initial movements");
        printf("\n");
        initialization = false;
      }
      fflush(stdout);
    }
    // ---------------------------------------------------------------
    // fine-tune scale so that the acceptance rate will be reasonable
    // ---------------------------------------------------------------
    if (iteration % 100 == 0) {
      const double acceptance_rate =
        (double)(workingmoves-oldworkingmoves)/(totalmoves-oldtotalmoves);
      oldworkingmoves = workingmoves;
      oldtotalmoves = totalmoves;
      if (acceptance_rate < 0.5) {
        scale /= 1.01;
        theta_scale /= 1.01;
      } else if (acceptance_rate > 0.7) {
        scale *= 1.01;
        theta_scale *= 1.01;
      }
      if (iteration % 1000 == 0)
        printf("scale: %g, acc: %g\n\n", scale, acceptance_rate);
    }
    iteration ++;
  }
  printf("---------------------------------------------------------------\n");
  printf("Initialization complete!\n");
  printf("---------------------------------------------------------------\n");

  // ---------------------------------------------------------------------------
  // MAIN PROGRAM LOOP
  // ---------------------------------------------------------------------------
  clock_t output_period = CLOCKS_PER_SEC; // start at outputting every minute fixme
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*60; // top out at one hour interval
  clock_t last_output = clock(); // when we last output data

  int frame = 0;
  totalmoves = 0;
  workingmoves = 0;
  for(long iteration=1; iteration<=iterations; iteration++) {
    // ---------------------------------------------------------------
    // Move each polyhedron once
    // ---------------------------------------------------------------
    totalmoves += N;
    workingmoves += move_all(polyhedra, N, periodic, walls, real_walls,
                             neighborR, scale, theta_scale, max_neighbors);
    // ---------------------------------------------------------------
    // Add data to historams
    // ---------------------------------------------------------------
    for(int i=0; i<N; i++) {
      // Density histogram:
      const int x_i = floor(polyhedra[i].pos[0]/de_density);
      const int y_i = floor(polyhedra[i].pos[1]/de_density);
      const int z_i = floor(polyhedra[i].pos[2]/de_density);
      density_histogram[x_i] ++;
      density_histogram[int(round(len[0]/de_density)) + y_i] ++;
      density_histogram[int(round((len[0] + len[1])/de_density)) + z_i] ++;
    }
    // ---------------------------------------------------------------
    // Save to file
    // ---------------------------------------------------------------
    const clock_t now = clock();
    if ((now > last_output + output_period) || iteration==iterations) {
      last_output = now;
      assert(last_output);
      if (output_period < max_output_period/2)
        output_period *= 2;
        else if (output_period < max_output_period)
          output_period = max_output_period;
      const double secs_done = double(now)/CLOCKS_PER_SEC;
      const int seconds = int(secs_done) % 60;
      const int minutes = int(secs_done / 60) % 60;
      const int hours = int(secs_done / 3600) % 24;
      const int days = int(secs_done / 86400);
      printf("Saving data after %i days, %02i:%02i:%02i, %li iterations complete.\n",
             days, hours, minutes, seconds, iteration);
      fflush(stdout);

      char *headerinfo = new char[4096];
      sprintf(headerinfo, "\
# period: (%5.2f, %5.2f, %5.2f), walls: (%5.2f, %5.2f, %5.2f), de_density: %g\n\
# seed: %li, R: %f, scale: %g, theta_scale: %g, real_walls: %i\n\
# iteration: %li, workingmoves: %li, totalmoves: %li, acceptance rate: %g\n",
              periodic[0], periodic[1], periodic[2],
              walls[0], walls[1], walls[2], de_density, seed, R,
              scale, theta_scale, real_walls, iteration, workingmoves,
              totalmoves, double(workingmoves)/totalmoves);

      // Saving density in each of the x, y, z dimensions
      char *density_filename = new char[1024];
      sprintf(density_filename, "%s/%s-density-%s-%i.dat", dir, filename, shape->name, N);
      FILE *densityout = fopen((const char *)density_filename, "w");
      delete[] density_filename;
      fprintf(densityout, "%s", headerinfo);
      fprintf(densityout, "\n#e       xdensity   ydensity   zdensity   histograms in x,y,z order\n");
      const int xbins = round(len[0]/de_density);
      const int ybins = round(len[1]/de_density);
      const int zbins = round(len[2]/de_density);
      int maxbin = max(max(xbins, ybins), zbins);
      for(int e_i = 0; e_i < maxbin; e_i ++) {
        const double e = (e_i + 0.5)*de_density;
        const double xshell_volume = len[1]*len[2]*de_density;
        const double yshell_volume = len[0]*len[2]*de_density;
        const double zshell_volume = len[0]*len[1]*de_density;
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
      delete[] headerinfo;
    }
    // ---------------------------------------------------------------
    // Save locations of vertices if desired
    // ---------------------------------------------------------------
    if (vertex_period > 0 && iteration % vertex_period == 0) {
      printf("Saving vertex locations. Frame: %i. Iteration: %li.\n", frame, iteration);
      char *vertices_fname = new char[1024];
      sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape->name, N, frame);
      save_locations(polyhedra, N, vertices_fname, iteration);
      delete[] vertices_fname;
      frame ++;
    }
  }
  // ---------------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // ---------------------------------------------------------------------------
  print_bad(polyhedra, N, periodic);

  for(int i=0; i<N; i++)
    delete[] polyhedra[i].neighbors;
  delete[] polyhedra;
  delete[] density_histogram;
  return 0;
}
// -----------------------------------------------------------------------------
// END OF MAIN
// -----------------------------------------------------------------------------






// Uses the seperating axis theorem. Not the fastest algorithm, but simple to implement.
// Theorem: Two convex objects do not overlap iff there exists a line onto which their
// 1d projections do not overlap.
// In three dimensions, if such a line exists, then the normal line to one of the
// faces of one of the shapes will be such a line.
inline bool overlap(const polyhedron &a, const polyhedron &b, double periodic[3]) {
  const vector3d ab = periodic_diff(a.pos, b.pos, periodic);
  if (ab.normsquared() > (a.R + b.R)*(a.R + b.R))
    return false;
  // construct axes from a
  // project a and b to each axis
  for (int i=0; i<a.mypoly->nfaces; i++) {
    const vector3d axis = a.rot.rotate_vector(a.mypoly->faces[i]);
    double projection = axis.dot(a.rot.rotate_vector((a.mypoly->vertices[0])*a.R));
    double mina = projection, maxa = projection;
    for (int j=1; j<a.mypoly->nvertices; j++) {
      projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*a.R));
      if (projection < mina) mina = projection;
      else if (projection > maxa) maxa = projection;
    }
    projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[0]*b.R) + ab);
    double minb = projection, maxb = projection;
    for (int j=1; j<b.mypoly->nvertices; j++) {
      projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[j]*b.R) + ab);
      if (projection < minb) minb = projection;
      else if (projection > maxb) maxb = projection;
    }
    if (mina > maxb || minb > maxa) {
      return false;
    }
  }
  // construct axes from b
  // project a and b to each axis
  for (int i=0; i<b.mypoly->nfaces; i++) {
    const vector3d axis = b.rot.rotate_vector(b.mypoly->faces[i]);
    double projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[0]*a.R));
    double mina = projection, maxa = projection;
    for (int j=1; j<a.mypoly->nvertices; j++) {
      projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*a.R));
      if (projection < mina) mina = projection;
      else if (projection > maxa) maxa = projection;
    }
    projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[0]*b.R) + ab);
    double minb = projection, maxb = projection;
    for (int j=1; j<b.mypoly->nvertices; j++) {
      projection = axis.dot(b.rot.rotate_vector(b.mypoly->vertices[j]*b.R) + ab);
      if (projection < minb) minb = projection;
      else if (projection > maxb) maxb = projection;
    }
    if (mina > maxb || minb > maxa) {
      return false;
    }
  }
  return true;
}

inline int overlaps_with_any(const polyhedron &a, const polyhedron *bs, double periodic[3], bool count, double dr) {
  // construct axes from a and a's projection onto them
  vector3d aaxes[a.mypoly->nfaces];
  double amins[a.mypoly->nfaces], amaxes[a.mypoly->nfaces];
  for (int i=0; i<a.mypoly->nfaces; i++) {
    aaxes[i] = a.rot.rotate_vector(a.mypoly->faces[i]);
    double projection = aaxes[i].dot(a.rot.rotate_vector(a.mypoly->vertices[0]*a.R));
    amins[i] = projection, amaxes[i] = projection;
    for (int j=1; j< a.mypoly->nvertices; j++) {
      projection = aaxes[i].dot(a.rot.rotate_vector(a.mypoly->vertices[j]*a.R));
      amins[i] = min(projection, amins[i]);
      amaxes[i] = max(projection, amaxes[i]);
    }
  }
  int num_overlaps = 0;
  for (int l=0; l<a.num_neighbors; l++) {
    const int k = a.neighbors[l];
    const vector3d ab = periodic_diff(a.pos, bs[k].pos, periodic);
    if (ab.normsquared() < (a.R + bs[k].R)*(a.R + bs[k].R)) {
      bool overlap = true; // assume overlap until we prove otherwise or fail to.
      // check projection of b against a's axes
      for (int i=0; i<a.mypoly->nfaces; i++) {
        double projection = aaxes[i].dot
          (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[0]*bs[k].R) + ab);
        double bmin = projection, bmax = projection;
        for (int j=1; j<bs[k].mypoly->nvertices; j++) {
          projection = aaxes[i].dot
            (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[j]*bs[k].R) + ab);
          bmin = min(projection, bmin);
          bmax = max(projection, bmax);
        }
        if (amins[i] > bmax || bmin > amaxes[i]) {
          overlap = false;
          i = a.mypoly->nfaces; // no overlap, move on to next
        }
      }
      if (overlap) { // still need to check against b's axes
        for (int i=0; i<bs[k].mypoly->nfaces; i++) {
          const vector3d axis = bs[k].rot.rotate_vector(bs[k].mypoly->faces[i]);
          double projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[0]*a.R));
          double amin = projection, amax = projection;
          for (int j=1; j<a.mypoly->nvertices; j++) {
            projection = axis.dot(a.rot.rotate_vector(a.mypoly->vertices[j]*a.R));
            amin = min(projection, amin);
            amax = max(projection, amax);
          }
          projection = axis.dot
            (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[0]*bs[k].R) + ab);
          double bmin = projection, bmax = projection;
          for (int j=1; j<bs[k].mypoly->nvertices; j++) {
            projection = axis.dot
              (bs[k].rot.rotate_vector(bs[k].mypoly->vertices[j]*bs[k].R) + ab);
            bmin = min(projection, bmin);
            bmax = max(projection, bmax);

          }
          if (amin > bmax || bmin > amax) {
            overlap = false;
            i = bs[k].mypoly->nfaces; //no overlap, move on to next
          }
        }
      }
      if (overlap) {
        if(!count)
          return 1;
        num_overlaps ++;
      }
    }
  }
  return num_overlaps;
}

inline bool in_cell(const polyhedron &p, double walls[3], bool real_walls) {
  if(real_walls) {
    for (int i=0; i<3; i++) {
      if (walls[i] > 0) {
        if (p.pos[i] - p.R > 0.0 && p.pos[i] + p.R < walls[i]) {
          continue;
        }
        double coord = (p.rot.rotate_vector(p.mypoly->vertices[0]*p.R) + p.pos)[i];
        double pmin = coord, pmax = coord;
        for (int j=1; j<p.mypoly->nvertices; j++) {
          coord = (p.rot.rotate_vector(p.mypoly->vertices[j]*p.R) + p.pos)[i];
          pmin = min(coord, pmin);
          pmax = max(coord, pmax);
        }
        if (pmin < 0.0 || pmax > walls[i])
          return false;
      }
    }
  }
  return true;
}

inline vector3d periodic_diff(const vector3d &a, const vector3d &b, double periodic[3]) {
  vector3d v = b - a;
  for (int i=0; i<3; i++) {
    if (periodic[i] > 0) {
      while (v[i] > periodic[i]/2.0)
        v[i] -= periodic[i];
      while (v[i] < -periodic[i]/2.0)
        v[i] += periodic[i];
    }
  }
  return v;
}

inline vector3d move(const vector3d &v, double scale, const double len[3]) {
  const vector3d newv = v + vector3d::ran(scale);
  return fix_periodic(newv, len);
}

inline vector3d fix_periodic(vector3d newv, const double len[3]) {
  for (int i=0; i<3; i++) {
    while (newv[i] > len[i])
      newv[i] -= len[i];
    while (newv[i] < 0.0)
      newv[i] += len[i];
    }
  return newv;
}

inline void update_neighbors(polyhedron &a, int n, const polyhedron *bs, int N, double neighborR, double periodic[3]) {
  a.num_neighbors = 0;
  for(int i=0; i<N; i++) {
    if ((i!=n) && (periodic_diff(a.pos, bs[i].neighbor_center, periodic).normsquared() <
         sqr(a.R + bs[i].R + neighborR))) {
      a.neighbors[a.num_neighbors] = i;
      a.num_neighbors ++;
    }
  }
}

inline void inform_neighbors(const polyhedron &newp, const polyhedron &oldp, int n, polyhedron *p, int N) {
  if (debug && false) {
    printf("informing: %i\n", n);
  }
  int new_index = 0, old_index = 0;

  // while (1) {
  //   if (new_index >= newlen) {
  //     handle olds;
  //     return;
  //   }
  //   if (old_index >= oldlen) {
  //     handle news;
  //     return;
  //   }
  //   if (oldindex >= oldlen || new[new_index] == old[old_index]) {
  //   } else if (new[index] < ) {
  //     handle new;
  //     new_index++;
  //   } else {
  //     handle old;
  //     old_index++;
  //   }
  // }
  bool checknew = newp.num_neighbors > 0;
  bool checkold = oldp.num_neighbors > 0;
  while (checknew || checkold) {
    if (debug && false) {
      if (checknew)
        printf("new[%i]: %i\t", new_index, newp.neighbors[new_index]);
      else
        printf("           \t");
      if (checkold)
        printf("old[%i]: %i", old_index, oldp.neighbors[old_index]);
      printf("\n");
    }
    if (checkold &&
        (!checknew || oldp.neighbors[old_index] < newp.neighbors[new_index])) {
      //printf("old[%i]: %i\n", old_index, oldp.neighbors[old_index]);
      // We had a neighbor that we moved away from, so it's time to remove ourselves
      // from his table and never speak to him again.
      const int i = oldp.neighbors[old_index];
      int pindex = p[i].num_neighbors - 1;
      int temp = p[i].neighbors[pindex];
      while (temp != n) {
        pindex --;
        const int temp2 = temp;
        temp = p[i].neighbors[pindex];
        p[i].neighbors[pindex] = temp2;
      }
      p[i].num_neighbors --;
      old_index ++;
    } else if (checknew &&
               (!checkold || newp.neighbors[new_index] < oldp.neighbors[old_index])) {
      // We have a new neighbor! Time to insert ourselves into his list
      // so we can be invited over for dinner.
      // We only increment the new index here.
      const int i = newp.neighbors[new_index];
      int pindex = p[i].num_neighbors - 1;
      int temp = p[i].neighbors[pindex];
      p[i].num_neighbors ++;
      while  (temp > n) {
        p[i].neighbors[pindex + 1] = temp;
        pindex --;
        temp = p[i].neighbors[pindex];
      }
      p[i].neighbors[pindex+1] = n;
      new_index ++;
    } else {// (newp.neighbors[new_index] == oldp.neighbors[old_index])
      // This new neighbor was also an old neighbor, so we don't have to do anything
      // but increment.
      new_index ++;
      old_index ++;
    }
    checknew = new_index < newp.num_neighbors;
    checkold = old_index < oldp.num_neighbors;
    if (debug) {
      if(new_index > newp.num_neighbors || old_index > oldp.num_neighbors) {
        printf("We have a problem. Inform neighbors not working as intended.\n\tNew index: %i, num new neighbors: %i, old index: %i, old num neighbors: %i.\n", new_index, newp.num_neighbors, old_index, oldp.num_neighbors);
        testcase = true;
      }
    }
  }
}

inline int move_all(polyhedron *p, int N, double periodic[3], double walls[3], bool real_walls, double neighborR, double scale, double theta_scale, int max_neighbors) {
  int workingmoves = 0;
  for(int i=0; i<N; i++) {
    if (debug && testcase)
      return -1;
    polyhedron temp;
    const double len[3] = {periodic[0]+walls[0], periodic[1]+walls[1], periodic[2]+walls[2]};
    temp.pos = move(p[i].pos, scale, len);
    temp.rot = rotation::ran(theta_scale)*p[i].rot;
    temp.R = p[i].R;
    temp.mypoly = p[i].mypoly;
    temp.neighbors = p[i].neighbors;
    temp.num_neighbors = p[i].num_neighbors;
    temp.neighbor_center = p[i].neighbor_center;
    if (in_cell(temp, walls, real_walls)) {
      bool overlaps = overlaps_with_any(temp, p, periodic);
      if (!overlaps) {
        const bool get_new_neighbors =
          (periodic_diff(temp.pos, p[i].neighbor_center, periodic).normsquared() >
           sqr(neighborR/2.0));
        if (get_new_neighbors) {
          // If we've moved too far, then the overlap test may have given a false
          // negative. So we'll find our new neighbors, and check against them.
          // If we still don't overlap, then we'll have to update the tables
          // of our neighbors that have changed.
          temp.neighbors = new int[max_neighbors];
          update_neighbors(temp, i, p, N, neighborR, periodic);
          if (debug) check_neighbor_symmetry(p, N);
          // However, for this check (and this check only), we don't need to
          // look at all of our neighbors, only our new ones.
          // fixme: do this!
          //int *new_neighbors = new int[max_neighbors];

          overlaps = overlaps_with_any(temp, p, periodic);
          if (!overlaps) {
            // Okay, we've checked twice, just like Santa Clause, so we're definitely
            // keeping this move and need to tell our neighbors where we are now.
            temp.neighbor_center = temp.pos;
            if (debug && testcase) print_one(p[i], -i, p, N, periodic);
            inform_neighbors(temp, p[i], i, p, N);
            delete[] p[i].neighbors;
          } else {
            // We don't have to move! We can throw away the new table as we don't
            // have any new neighbors.
            delete[] temp.neighbors;
          }
        }
        if (!overlaps) {
          p[i] = temp;
          workingmoves ++;
        }
      }
      // -----------------------------------------------------------------------
      // debug tests
      // -----------------------------------------------------------------------
      if (debug) {
        bool test_overlap = false;
        for(int j=0; j<N; j++) {
          if (j != i && overlap(temp, p[j], periodic)) {
            test_overlap = true;
            break;
          }
        }
        if (test_overlap != overlaps) {
          printf("\n\nError error! We do not have agreement on whether or not they overlap!\ndebug says: %i but normal test says: %i for shape: %i.\n\n", test_overlap, overlaps, i);
          p[i] = temp;
          print_all(p, N, periodic);
          // char *vertices_fname = new char[1024];
          // sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape->name, N, -1);
          // save_locations(p, vertices_fname);
          // delete[] vertices_fname;
          testcase = true;
          return -1;
        }
      }
      // -----------------------------------------------------------------------
      // done with debug tests
      // -----------------------------------------------------------------------
    }
  }
  return workingmoves;
}




inline void print_all(const polyhedron *p, int N, double periodic[3]) {
  if (debug) {
    for (int i=0; i<N; i++) {
      char *pos = new char[1024];
      char *rot = new char[1024];
      p[i].pos.tostr(pos);
      p[i].rot.tostr(rot);
      printf("%4i: %s, R: %4.2f, %i neighbors: ", i, p[i].mypoly->name, p[i].R, p[i].num_neighbors);
      for(int j=0; j<min(10, p[i].num_neighbors); j++)
        printf("%i ", p[i].neighbors[j]);
      if (p[i].num_neighbors > 10)
        printf("...");
      printf("\n      pos:          %s\n      rot: %s\n", pos, rot);
    }
    printf("\n");
    fflush(stdout);
  }
}

inline void print_one(const polyhedron &a, int id, const polyhedron *p, int N, double periodic[3]) {
  if (debug) {
    char *pos = new char[1024];
    char *rot = new char[1024];
    a.pos.tostr(pos);
    a.rot.tostr(rot);
    printf("%4i: %s, R: %4.2f, %i neighbors: ", id, a.mypoly->name, a.R, a.num_neighbors);
    for(int j=0; j<min(10, a.num_neighbors); j++)
      printf("%i ", a.neighbors[j]);
    if (a.num_neighbors > 10)
      printf("...");
    printf("\n      pos:          %s\n      rot: %s\n", pos, rot);
    for (int j=0; j<N; j++) {
      if (j != id && overlap(a, p[j], periodic)) {
        p[j].pos.tostr(pos);
        p[j].rot.tostr(rot);
        printf("\t  Overlaps with %i", j);
        printf(": %s   %s\n", pos, rot);
      }
    }
  }
  printf("\n");
  fflush(stdout);
}

inline void print_bad(const polyhedron *p, int N, double periodic[3]) {
  if (debug) {
    for (int i=0; i<N; i++) {
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
        char *rot = new char[1024];
        p[i].pos.tostr(pos);
        p[i].rot.tostr(rot);
        printf("%s %4i: %s   %s R: %4.2f\n", p[i].mypoly->name, i, pos, rot, p[i].R);
        if (!incell)
          printf("\t  Outside cell!\n");
        for (int j=0; j<N; j++) {
          if (j != i && overlap(p[i], p[j], periodic)) {
            p[j].pos.tostr(pos);
            p[j].rot.tostr(rot);
            printf("\t  Overlaps with %i", j);
            printf(": %s   %s\n", pos, rot);
          }
        }
        delete[] pos;
        delete[] rot;
      }
    }
  }
  fflush(stdout);
}

inline void check_neighbor_symmetry(const polyhedron *p, int N) {
  if (debug) {
    for(int i=0; i<N; i++) {
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

void save_locations(const polyhedron *p, int N, const char *fname, long iteration) {
  FILE *out = fopen((const char *)fname, "w");
  fprintf(out, "# iteration: %li.\n", iteration);
  for(int i=0; i<N; i++) {
    fprintf(out, "%6.2f %6.2f %6.2f ", p[i].pos[0], p[i].pos[1], p[i].pos[2]);
    for(int j=0; j<p[i].mypoly->nvertices; j++) {
      const vector3d vertex =
        p[i].rot.rotate_vector(p[i].mypoly->vertices[j]*p[i].R) + p[i].pos;
      fprintf(out, "%6.2f %6.2f %6.2f ", vertex[0], vertex[1], vertex[2]);
    }
    fprintf(out, "\n");
  }
  fclose(out);
}

poly_shape::poly_shape(const char *set_name) {
  if (strcmp(set_name, "cube") == 0) {
    nvertices = 8;
    nfaces = 3;
    volume = 1.0;
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    name = new char[4];
    sprintf(name, "cube");

    const double v_cube = 1.0/sqrt(3.0);
    vertices[0] = vector3d( v_cube,  v_cube,  v_cube);
    vertices[1] = vector3d(-v_cube,  v_cube,  v_cube);
    vertices[2] = vector3d(-v_cube, -v_cube,  v_cube);
    vertices[3] = vector3d( v_cube, -v_cube,  v_cube);
    vertices[4] = vector3d( v_cube,  v_cube, -v_cube);
    vertices[5] = vector3d(-v_cube,  v_cube, -v_cube);
    vertices[6] = vector3d(-v_cube, -v_cube, -v_cube);
    vertices[7] = vector3d( v_cube, -v_cube, -v_cube);

    faces[0] = vector3d(1, 0, 0);
    faces[1] = vector3d(0, 1, 0);
    faces[2] = vector3d(0, 0, 1);
  }
  else if (strcmp(set_name, "tetrahedron") == 0) {
    nvertices = 4;
    nfaces = 4;
    volume = 1.0/3.0;
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    name = new char[11];
    sprintf(name, "tetrahedron");

    vertices[0] = vector3d( sqrt(2.0/3.0),    -sqrt(2.0)/3.0, -1.0/3.0);
    vertices[1] = vector3d(-sqrt(2.0/3.0),    -sqrt(2.0)/3.0, -1.0/3.0);
    vertices[2] = vector3d(             0, 2.0*sqrt(2.0)/3.0, -1.0/3.0);
    vertices[3] = vector3d(             0,                 0,      1.0);

    faces[0] = vector3d(         0, -sqrt(2.0),           0.5);
    faces[1] = vector3d( sqrt(3.0),        1.0, 1.0/sqrt(2.0));
    faces[2] = vector3d(-sqrt(3.0),        1.0, 1.0/sqrt(2.0));
    faces[3] = vector3d(         0,          0,             1);
  }
  else if (strcmp(set_name, "truncated_tetrahedron") == 0) {
    nvertices = 18;
    nfaces = 4;
    volume = 0; //fixme
    vertices = new vector3d[nvertices];
    faces = new vector3d[nfaces];
    name = new char[21];
    sprintf(name, "truncated_tetrahedron");

    faces[0] = vector3d(         0, -sqrt(2.0),           0.5);
    faces[1] = vector3d( sqrt(3.0),        1.0, 1.0/sqrt(2.0));
    faces[2] = vector3d(-sqrt(3.0),        1.0, 1.0/sqrt(2.0));
    faces[3] = vector3d(         0,          0,             1);
  }
  else {
    nvertices = 0;
    nfaces = 0;
    volume = 0;
    vertices = NULL;
    faces = NULL;
    name = new char[13];
    sprintf(name, "invalid shape");
  }
}
