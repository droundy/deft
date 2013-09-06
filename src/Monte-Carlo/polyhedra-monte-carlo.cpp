#include <stdio.h>
#include <time.h>
#include "MersenneTwister.h"
#include "vector3d.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utilities.h"
#include "handymath.h"

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

  poly_shape(int num_verts, int num_faces, const char *set_name) {
    nvertices = num_verts;
    nfaces = num_faces;
    vertices = new vector3d[nvertices]();
    faces = new vector3d[nfaces]();
    name = new char[1024]();
    sprintf(name, "%s", set_name);
  }
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
  poly_shape *mypoly;
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

// We only need to check overlaps with neighbors, but we need to update
// neighbor tables as soon as a polyhedron moves more than a distance of
// neighborR/2.0.
const double neighborR = 0.5;

// -----------------------------------------------------------------------------
// Global "Constants" -- set at runtime then unchanged
// -----------------------------------------------------------------------------
long unsigned int seed = 0; // for random number generator
int iterations;
int N;
bool periodic[3] = {false, false, false};
bool walls[3] = {false, false, false};
bool real_walls = true; // true means shapes can't overlap wall at all,
                        // false means just centers can't overlap wall
bool save_vertices = false;
int vertex_save_per_iters;
double len[3] = {20, 20, 20};
double de_density = 0.1;
double R = 1.0;
double dr = 0.01;
double scale = 0.005;
double theta_scale = 0.005;


// The most neighbors that any polyhedron could have.
// Not an exact value, but an over approximation based on volumes.
// Used for length of neighbors arrays
int max_neighbors;

// All of the shapes to use:
poly_shape cube(8, 3, "cube");
poly_shape tetrahedron(4, 4, "tetrahedron");
poly_shape truncated_tetrahedron(18, 4, "truncated tetrahedron");

// -----------------------------------------------------------------------------
// Global variables -- try to avoid making any more of these
// -----------------------------------------------------------------------------
long neighbor_updates = 0, neighbor_informs = 0;
long totalmoves = 0, workingmoves = 0;


// for debugging purposes
bool testcase = false;

// -----------------------------------------------------------------------------
// Functions
// -----------------------------------------------------------------------------

// Sets up the vertices, faces, etc. for each polyhedron. Only run once.
void define_shapes();

// Find's the neighbors of a by comparing a's position to the center of
// everyone else's neighborsphere, where n is the index of a in bs.
inline void update_neighbors(polyhedron &a, const polyhedron *bs, int n);

// Removes n from the neighbor table of anyone neighboring oldp but not newp.
// Adds n to the neighbor table of anyone neighboring newp but not oldp.
// Keeps the neighbor tables sorted.
inline void inform_neighbors(const polyhedron &newp, const polyhedron &oldp, int n, polyhedron *p);

// Attempts to move each polyhedron once. Updates the global variables
// totalmoves and workingmoves.
inline void move_all(polyhedron *p);

// Check whether two polyhedra overlap
inline bool overlap(const polyhedron &a, const polyhedron &b);

// Check whether a polyhedron overlaps with any in the array, except for
// Useful for post-move testing with a temporary polyhedron.
// If count is true, it will return the total number of overlaps, otherwise
// it returns 1 if there is at least one overlap, 0 if there are none.
inline int overlaps_with_any(const polyhedron &a, const polyhedron *bs, bool count=false, double dr=0);

// Check if polyhedron is outside its allowed locations
inline bool in_cell(const polyhedron &p);

// Return the vector pointing from a to b, accounting for periodic
// boundaries
inline vector3d periodic_diff(const vector3d &a, const vector3d  &b);

// Move v in a random direction by a distance determined by a gaussian distribution
inline vector3d move(const vector3d &v, double scale);

// If v is outside the cell, and there are periodic boundary condition(s), it is
// moved (appropriately) into the cell
inline vector3d fix_periodic(vector3d newv);

// States how long it's been since last took call.
static void took(const char *name);

// Saves the vertices of all polyhedra to a file.
void save_locations(const polyhedron *p, const char *fname);



// The following functions only do anything if debug is true:

// Prints the locations, rotations, radii, and shapes of every polyhedron
// As well as if any overlap or are outside the cell.
inline void print_all(const polyhedron *p);

// Same as print_all, but only prints information for one polyhedron,
// and does not do overlap tests.
inline void print_one(const polyhedron &p, int id);

// Only print those shapes that overlap or are outside the cell
// Also prints those they overlap with
inline void print_bad(const polyhedron *p);

// Checks to make sure that every polyhedron is his neighbor's neighbor.
inline void check_neighbor_symmetry(const polyhedron *p);



int main(int argc, const char *argv[]) {
  // Setup all the basic types of polyhedra
  define_shapes();
  // create a generic poly_shape that each polyhedron will point to
  // (assuming we don't have a mixture)
  poly_shape *shape = &cube;
  // ---------------------------------------------------------------------------
  // Set values from parameters
  // ---------------------------------------------------------------------------
  const int minparams = 5;
  if (argc < minparams) {
    printf("\nUse: %s N iterations directory filename\n", argv[0]);
    printf("Additional optional parameters:\n\t\tperiodx, periody, periodz, wallx, wally, wallz,\n");
    printf("\t\tdimensions {LENX LENY LENZ}, de_density {VAL}, seed {VAL}, shape {POLY},\n");
    printf("\t\tR {VAL}, scale {VAL}, theta_scale {VAL}, fake_walls, save_vertices {ITERS}\n");
    printf("Anything in brackets means to put values there.\n");
    printf("Available shapes currently include: cube and tetrahedron.\n");
    return 11;
  }
  N = atoi(argv[1]);
  const long iterations = atol(argv[2]);
  const char *dir = argv[3];
  const char *filename = argv[4];
  printf("----------------------------------------------------------------------\n");
  printf("Running %s with the following parameters:\n", argv[0]);
  printf("Setting number of polyhedra to %i\n", N);
  printf("Running for %4.2e iterations\n", double(iterations));
  printf("Saving data to: %s/%s\n", dir, filename);
  for (int i=minparams; i<argc; i++) {
    if (strcmp(argv[i], "dimensions") == 0) {
      len[0] = atof(argv[i+1]);
      len[1] = atof(argv[i+2]);
      len[2] = atof(argv[i+3]);
      printf("Setting cell dimensions to (%g, %g, %g).\n", len[0], len[1], len[2]);
      i += 3;
    } else if (strcmp(argv[i], "periodx") == 0) {
      printf("Turning on periodic in the x dimension.\n");
      periodic[0] = true;
    } else if (strcmp(argv[i], "periody") == 0) {
      printf("Turning on periodic in the y dimension.\n");
      periodic[1] = true;
    } else if (strcmp(argv[i], "periodz") == 0) {
      printf("Turning on periodic in the z dimension.\n");
      periodic[2] = true;
    } else if (strcmp(argv[i], "wallx") == 0) {
      printf("Turning on walls in the x dimension.\n");
      walls[0] = true;
    } else if (strcmp(argv[i], "wally") == 0) {
      printf("Turning on walls in the y dimension.\n");
      walls[1] = true;
    } else if (strcmp(argv[i], "wallz") == 0) {
      printf("Turning on walls in the z dimension.\n");
      walls[2] = true;
    } else if (strcmp(argv[i], "de_density") == 0) {
      de_density = atof(argv[i+1]);
      printf("Setting density resolution to %g.\n", de_density);
      i += 1;
    } else if (strcmp(argv[i], "R") == 0) {
      R = atof(argv[i+1]);
      printf("Setting radius of the circumscribing sphere to %g.\n", R);
      i += 1;
    } else if (strcmp(argv[i], "dr") == 0) {
      dr = atof(argv[i+1]);
      printf("Setting radius of the circumscribing sphere to %g.\n", dr);
      i += 1;
    } else if (strcmp(argv[i], "scale") == 0) {
      scale = atof(argv[i+1]);
      printf("Setting movement scale to %g.\n", scale);
      i += 1;
    } else if (strcmp(argv[i], "theta_scale") == 0) {
      theta_scale = atof(argv[i+1]);
      printf("Setting rotation scale to %g.\n", theta_scale);
      i += 1;
    } else if (strcmp(argv[i], "seed") == 0) {
      seed = atol(argv[i+1]);
      printf("Setting random seed to %lu.\n", seed);
      i += 1;
    } else if (strcmp(argv[i], "fake_walls") == 0) {
      real_walls = false;
      printf("Using fake walls. Intersections will only be determined based on centers of polyhedra.\n");
    } else if (strcmp(argv[i], "save_vertices") == 0) {
      save_vertices = true;
      vertex_save_per_iters = atoi(argv[i+1]);
      printf("Will save the coordinates of the vertices of all shapes every %i iterations.\n", vertex_save_per_iters);
      i += 1;
    } else if (strcmp(argv[i], "shape") == 0) {
      if (strcmp(argv[i+1], "cube") == 0) {
        shape = &cube;
        printf("Using cubes.\n");
      } else if (strcmp(argv[i+1], "tetrahedron") == 0) {
        shape = &tetrahedron;
        printf("Using cubes.\n");
      } else {
        printf("Invalid shape: %s.\n", argv[i+1]);
        return 13;
      }
      i += 1;
    } else {
      printf("Invalid parameter: %s.\n", argv[i]);
      return 15;
    }
  }
  printf("Neighbor radius is %g.\n", neighborR);
  if (totime) printf("Timing information will be displayed.\n");
  if (debug) printf("DEBUG MODE IS ENABLED!\n");
  else printf("Debug mode disabled.\n");
  printf("----------------------------------------------------------------------\n\n");

  // Create arrays
  const int density_bins = round((len[0] + len[1] + len[2])/de_density);
  long *density_histogram = new long[3*density_bins]();
  polyhedron *polyhedra = new polyhedron[N];

  // Initialize the random number generator with our seed so we get repreducable
  // results.

  ////////////////////////////////////////////////////////////////////////////////
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
  max_neighbors = min(N, neighbor_sphere_vol / shape->volume);

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
  printf("Setting up grid with polyhedra: (%i, %i, %i) and space: (%g, %g, %g).\n",
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

      polyhedra[i].pos[0] = (x)*xspace;
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
        (periodic_diff(polyhedra[i].pos, polyhedra[j].pos).normsquared() <
         uipow(polyhedra[i].R + polyhedra[j].R + neighborR, 2));
      if (is_neighbor) {
        const int index = polyhedra[i].num_neighbors;
        polyhedra[i].num_neighbors ++;
        polyhedra[i].neighbors[index] = j;
      }
    }
    if (debug)
      printf("Polyhedron %i has %i neighbors.\n", i, polyhedra[i].num_neighbors);
  }

  print_all(polyhedra);
  print_bad(polyhedra);
  printf("\n\n");
  // Make sure none are overlapping:
  for(int i=0; i<N; i++) {
    if (!in_cell(polyhedra[i])) {
      printf("Oops, this is embarassing. I seem to have placed some things outside our cell.\n");
      printf("You might want to look into that.\n");
      char *vertices_fname = new char[1024];
      sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape->name, N, -1);
      save_locations(polyhedra, vertices_fname);
      delete[] vertices_fname;
      return 17;
    }
    for(int j=i+1; j<N; j++) {
      if (overlap(polyhedra[i], polyhedra[j])) {
        printf("ERROR in intial placement. We have overlaps!!!\n");
        printf("AHHHHHH I DON'T KNOW WHAT TO DO!@!!!!1111\n");
        char *vertices_fname = new char[1024];
        sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape->name, N, -1);
        save_locations(polyhedra, vertices_fname);
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
  printf("// ---------------------------------------------------------------\n");
  printf("Beginning the initialization.\n");
  printf("// ---------------------------------------------------------------\n");
  long iteration = 1;
  long oldtotalmoves = 0, oldworkingmoves = 0;
  bool initialization = true;
  while (initialization) {
    if (debug) {
      print_bad(polyhedra);
    }
    // ---------------------------------------------------------------
    // Move each polyhedron once
    // ---------------------------------------------------------------
    move_all(polyhedra);

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
    iteration ++;
  }
  printf("// ---------------------------------------------------------------\n");
  printf("Initialization complete!\n");
  printf("// ---------------------------------------------------------------\n");

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
    move_all(polyhedra);
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
      sprintf(headerinfo, "# dim: (%5.2f, %5.2f, %5.2f), period: (%i, %i, %i), \
walls: (%i, %i, %i), de_density: %g, seed: %li\n# R: %f, scale: %g, theta_scale: %g, \
real_walls: %i\n# iteration: %li, workingmoves: %li, totalmoves: %li, acceptance \
rate: %g\n",
              len[0], len[1], len[2], periodic[0], periodic[1], periodic[2],
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
    if (save_vertices && iteration % vertex_save_per_iters == 0) {
      printf("Saving vertex locations. Frame: %i. Iteration: %li.\n", frame, iteration);
      char *vertices_fname = new char[1024];
      sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape->name, N, frame);
      save_locations(polyhedra, vertices_fname);
      delete[] vertices_fname;
      frame ++;
    }
  }
  // ---------------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // ---------------------------------------------------------------------------
  print_bad(polyhedra);

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
inline bool overlap(const polyhedron &a, const polyhedron &b) {
  const vector3d ab = periodic_diff(a.pos, b.pos);
  if (ab.normsquared() > (a.R + b.R)*(a.R + b.R))
    return false;
  // construct axes from a
  // project a and b to each axis
  for (int i=0; i<a.mypoly->nfaces; i++) {
    const vector3d axis = a.rot.rotate_vector(a.mypoly->faces[i]);
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

inline int overlaps_with_any(const polyhedron &a, const polyhedron *bs, bool count, double dr) {
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
    const vector3d ab = periodic_diff(a.pos, bs[k].pos);
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

inline bool in_cell(const polyhedron &p) {
  if(real_walls) {
    for (int i=0; i<3; i++) {
      if (walls[i]) {
        if (p.pos[i] - p.R > 0.0 && p.pos[i] + p.R < len[i]) {
          continue;
        }
        double coord = (p.rot.rotate_vector(p.mypoly->vertices[0]*p.R) + p.pos)[i];
        double pmin = coord, pmax = coord;
        for (int j=1; j<p.mypoly->nvertices; j++) {
          coord = (p.rot.rotate_vector(p.mypoly->vertices[j]*p.R) + p.pos)[i];
          pmin = min(coord, pmin);
          pmax = max(coord, pmax);
        }
        if (pmin < 0.0 || pmax > len[i])
          return false;
      }
    }
  }
  return true;
}

inline vector3d periodic_diff(const vector3d &a, const vector3d &b) {
  vector3d v = b - a;
  for (int i=0; i<3; i++) {
    if (periodic[i]) {
      while (v[i] > len[i]/2.0)
        v[i] -= len[i];
      while (v[i] < -len[i]/2.0)
        v[i] += len[i];
    }
  }
  return v;
}

inline vector3d move(const vector3d &v, double scale) {
  const vector3d newv = v + vector3d::ran(scale);
  return fix_periodic(newv);
}

inline vector3d fix_periodic(vector3d newv) {
  for (int i=0; i<3; i++) {
    if (periodic[i] || walls[i]) {
      while (newv[i] > len[i])
        newv[i] -= len[i];
      while (newv[i] < 0.0)
        newv[i] += len[i];
    }
  }
  return newv;
}

inline void update_neighbors(polyhedron &a, const polyhedron *bs, int n) {
  bool report = false;
  if (debug && (n == 0 || n == 6) && totalmoves > 6380 && false) {
    report = true;
    printf("updating: %i, r = %g\n", n, (R + R + neighborR));
  }
  a.num_neighbors = 0;
  for(int i=0; i<N; i++) {
    if ((i!=n) && (periodic_diff(a.pos, bs[i].neighbor_center).normsquared() <
         sqr(a.R + bs[i].R + neighborR))) {
      a.neighbors[a.num_neighbors] = i;
      a.num_neighbors ++;
    }
  }
  if (debug && report) {
    printf("%i has %i neighbors.\n", n, a.num_neighbors);
    for(int i=0; i<N; i++){
      char *pos = new char[1024];
      bs[i].neighbor_center.tostr(pos);
      printf("%i is %g from %i's neighborcenter at %s.\n", n, periodic_diff(a.pos, bs[i].neighbor_center).norm(), i, pos);
      delete[] pos;
    }
    print_one(a, n);
    printf("exiting update_neighbors.\n");
    if (totalmoves > 6385)
      testcase = true;
  }
}

inline void inform_neighbors(const polyhedron &newp, const polyhedron &oldp, int n, polyhedron *p) {
  if (debug && false) {
    printf("informing: %i\n", n);
    print_all(p);
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

inline void move_all(polyhedron *p) {
  for(int i=0; i<N; i++) {
    if (debug && testcase)
      return;
    totalmoves ++;
    polyhedron temp;
    temp.pos = move(p[i].pos, scale);
    temp.rot = rotation::ran(theta_scale)*p[i].rot;
    temp.R = p[i].R;
    temp.mypoly = p[i].mypoly;
    temp.neighbors = p[i].neighbors;
    temp.num_neighbors = p[i].num_neighbors;
    temp.neighbor_center = p[i].neighbor_center;
    if (in_cell(temp)) {
      bool overlaps = overlaps_with_any(temp, p);
      if (!overlaps) {
        const bool get_new_neighbors =
          (periodic_diff(temp.pos, p[i].neighbor_center).normsquared() >
           sqr(neighborR/2.0));
        if (get_new_neighbors) {
          // If we've moved too far, then the overlap test may have given a false
          // negative. So we'll find our new neighbors, and check against them.
          // If we still don't overlap, then we'll have to update the tables
          // of our neighbors that have changed.
          temp.neighbors = new int[max_neighbors];
          update_neighbors(temp, p, i);
          if (debug) check_neighbor_symmetry(p);
          neighbor_updates ++;
          // However, for this check (and this check only), we don't need to
          // look at all of our neighbors, only our new ones.
          // fixme: do this!
          //int *new_neighbors = new int[max_neighbors];

          overlaps = overlaps_with_any(temp, p);
          if (!overlaps) {
            // Okay, we've checked twice, just like Santa Clause, so we're definitely
            // keeping this move and need to tell our neighbors where we are now.
            temp.neighbor_center = temp.pos;
            if (debug && testcase) print_one(p[i], -i);
            inform_neighbors(temp, p[i], i, p);
            neighbor_informs ++;
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
          if (j != i && overlap(temp, p[j])) {
            test_overlap = true;
            break;
          }
        }
        if (test_overlap != overlaps) {
          printf("\n\nError error! We do not have agreement on whether or not they overlap!\ndebug says: %i but normal test says: %i for shape: %i.\n\n", test_overlap, overlaps, i);
          p[i] = temp;
          print_all(p);
          // char *vertices_fname = new char[1024];
          // sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape->name, N, -1);
          // save_locations(p, vertices_fname);
          // delete[] vertices_fname;
          testcase = true;
          return;
        }
      }
      // -----------------------------------------------------------------------
      // done with debug tests
      // -----------------------------------------------------------------------
    }
  }
}




inline void print_all(const polyhedron *p) {
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
      if (!in_cell(p[i]))
        printf("----  Outside cell!\n");
      for (int j=i+1; j<N; j++)
        if (overlap(p[i], p[j]))
          printf("****  Overlaps with %i!!!.\n", j);
    }
    printf("\n");
    fflush(stdout);
  }
}

inline void print_one(const polyhedron &p, int id) {
  if (debug) {
    char *pos = new char[1024];
    char *rot = new char[1024];
    p.pos.tostr(pos);
    p.rot.tostr(rot);
    printf("%4i: %s, R: %4.2f, %i neighbors: ", id, p.mypoly->name, p.R, p.num_neighbors);
    for(int j=0; j<min(10, p.num_neighbors); j++)
      printf("%i ", p.neighbors[j]);
    if (p.num_neighbors > 10)
      printf("...");
    printf("\n      pos:          %s\n      rot: %s\n", pos, rot);
    if (!in_cell(p))
      printf("----  Outside cell!\n");
  }
  printf("\n");
  fflush(stdout);
}

inline void print_bad(const polyhedron *p) {
  if (debug) {
    for (int i=0; i<N; i++) {
      bool incell = in_cell(p[i]), overlaps = false;
      for (int j=0; j<N; j++) {
        if (j != i && overlap(p[i], p[j])) {
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
          if (j != i && overlap(p[i], p[j])) {
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

inline void check_neighbor_symmetry(const polyhedron *p) {
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
          print_all(p);
          testcase = true;
        }
      }
    }
    if (testcase)
      print_all(p);
  }
}
static void took(const char *name) {
  assert(name); // so it'll count as being used...
  static clock_t last_time = clock();
  clock_t t = clock();
  double peak = peak_memory()/1024.0/1024;
  double seconds = (t-last_time)/double(CLOCKS_PER_SEC);
  if (seconds > 120) {
    printf("%s took %.0f minutes and %g M memory.\n", name, seconds/60, peak);
  } else {
    printf("%s took %g seconds and %g M memory.\n", name, seconds, peak);
  }
  fflush(stdout);
  last_time = t;
}

void save_locations(const polyhedron *p, const char *fname) {
  FILE *out = fopen((const char *)fname, "w");
  fprintf(out, "# saved every %i iterations, totalmoves: %li.\n", vertex_save_per_iters, totalmoves);
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

void define_shapes() {
  // Shapes, unrotated, at the origin, and circumscribed by a sphere of radius 1.0

  cube.volume = 1.0;
  const double v_cube = 1.0/sqrt(3.0);
  cube.vertices[0] = vector3d( v_cube,  v_cube,  v_cube);
  cube.vertices[1] = vector3d(-v_cube,  v_cube,  v_cube);
  cube.vertices[2] = vector3d(-v_cube, -v_cube,  v_cube);
  cube.vertices[3] = vector3d( v_cube, -v_cube,  v_cube);
  cube.vertices[4] = vector3d( v_cube,  v_cube, -v_cube);
  cube.vertices[5] = vector3d(-v_cube,  v_cube, -v_cube);
  cube.vertices[6] = vector3d(-v_cube, -v_cube, -v_cube);
  cube.vertices[7] = vector3d( v_cube, -v_cube, -v_cube);

  cube.faces[0] = vector3d(1, 0, 0);
  cube.faces[1] = vector3d(0, 1, 0);
  cube.faces[2] = vector3d(0, 0, 1);


  // Tetrahedron and truncated tetrahedron are oriented so that they have
  // a face parallel to the xy-plane and an edge parallel to the x-axis.
  tetrahedron.volume = 1.0/3.0;
  tetrahedron.vertices[0] = vector3d( sqrt(2.0/3.0),    -sqrt(2.0)/3.0, -1.0/3.0);
  tetrahedron.vertices[1] = vector3d(-sqrt(2.0/3.0),    -sqrt(2.0)/3.0, -1.0/3.0);
  tetrahedron.vertices[2] = vector3d(             0, 2.0*sqrt(2.0)/3.0, -1.0/3.0);
  tetrahedron.vertices[3] = vector3d(             0,                 0,      1.0);


  tetrahedron.faces[0] = vector3d(         0, -sqrt(2.0),           0.5);
  tetrahedron.faces[1] = vector3d( sqrt(3.0),        1.0, 1.0/sqrt(2.0));
  tetrahedron.faces[2] = vector3d(-sqrt(3.0),        1.0, 1.0/sqrt(2.0));
  tetrahedron.faces[3] = vector3d(         0,          0,             1);



  truncated_tetrahedron.faces[0] = vector3d(         0, -sqrt(2.0),           0.5);
  truncated_tetrahedron.faces[1] = vector3d( sqrt(3.0),        1.0, 1.0/sqrt(2.0));
  truncated_tetrahedron.faces[2] = vector3d(-sqrt(3.0),        1.0, 1.0/sqrt(2.0));
  truncated_tetrahedron.faces[3] = vector3d(         0,          0,             1);


}
