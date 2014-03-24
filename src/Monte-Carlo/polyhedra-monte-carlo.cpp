#include <stdio.h>
#include <time.h>
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <popt.h>
#include "handymath.h"
#include "vector3d.h"
#include "Monte-Carlo/polyhedra.h"


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

// Saves the vertices of all polyhedra to a file.
inline void save_locations(const polyhedron *p, int N, const char *fname,
                           const double len[3], const char *comment="");

// Generates all of the figures for the talk instead of performing a normal run
// returns 0 unless there's an error
inline int generate_talk_figs();

// The following functions only do anything if debug is true:

// Prints the locations, rotations, radii, and shapes of every polyhedron
// As well as if any overlap or are outside the cell.
inline void print_all(const polyhedron *p, int N, double periodic[3]);

// Same as print_all, but only prints information for one polyhedron,
// and does not do overlap tests.
inline void print_one(const polyhedron &a, int id, const polyhedron *p, int N,
                      double periodic[3]);

// Only print those shapes that overlap or are outside the cell
// Also prints those they overlap with
inline void print_bad(const polyhedron *p, int N, double periodic[3]);

// Checks to make sure that every polyhedron is his neighbor's neighbor.
inline void check_neighbor_symmetry(const polyhedron *p, int N);

int main(int argc, const char *argv[]) {
  took("Starting program");
  // -----------------------------------------------------------------------------
  // Define "Constants" -- set from arguments then unchanged
  // -----------------------------------------------------------------------------
  double periodic[3] = {0, 0, 0};
  double walls[3] = {0, 0, 0};
  bool fake_walls = false;
  unsigned long int seed = 0;
  int vertex_period = 0;
  char *shape_name = new char[1024];
  sprintf(shape_name, "truncated_tetrahedron");

  char *dir = new char[1024];
  sprintf(dir, "papers/polyhedra/figs/mc");
  char *filename = new char[1024];
  sprintf(filename, "[walls/periodic]-FF");
  char *structure = new char[1024];
  sprintf(structure, "ice");
  int N = 0;
  long iterations = 100000000000;
  long initialize_iterations = 500000;
  double acceptance_goal = .4;
  double R = 0;
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
    {"N", 'N', POPT_ARG_INT, &N, 0, "Number of polyhedra to simulate", "N"},
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
information on the number and type of polyhedra, as well as file extensions", "name"},
    {"dir", 'd', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &dir, 0,
     "Directory to save to", "dir"},
    {"shape", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &shape_name, 0,
     "Type of polyhedra to use. Can be one of [cube | tetrahedron | truncated_tetrahedron | cuboid]", "shape"},
    {"ratio", 'r', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &ratio, 0, "Ratio of sides of cuboid. For a cuboid of sides AxAxB, ratio = B/A.", "ratio"},
    {"R", 'R', POPT_ARG_DOUBLE, &R, 0,
     "Size of the sphere that circumscribes each polyhedron. Defaults to setting edge length to 1", "R"},
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
     "Standard deviation for translations of polyhedra", "scale"},
    {"theta_scale", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &theta_scale, 0,
     "Standard deviation for rotations of polyhedra", "theta"},
    {"seed", 's', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0,
     "Seed for the random number generator", "seed"},
    {"save_vertices", 'v', POPT_ARG_INT, &vertex_period, 0,
     "Periodically saves the vertex locations of the polyhedra, where period is a number \
of iterations", "period"},
    {"time", 't', POPT_ARG_INT, &totime, 0,
     "Timing information will be displayed", "interval"},
    {"acceptance_goal", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
     &acceptance_goal, 0, "Goal to set the acceptance rate", "goal"},
    {"talk", '\0', POPT_ARG_NONE, &talk, 0,
     "Generate figures for talk. Will not perform a normal run."},
    {"structure", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &structure, 0,
    "Structure to use for truncated tetrahedra. Does nothing for other shapes. \
Can be one of [ice | arsenic]"},
    POPT_AUTOHELP
    POPT_TABLEEND
  };
  optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTION...]\nNumber of polyhedra and \
periodicity/dimensions are not set by default and are required.\n");

  int c = 0;
  // go through arguments, set them based on optionsTable
  while((c = poptGetNextOpt(optCon)) >= 0);
  if (c < -1) {
    fprintf(stderr, "\n%s: %s\n", poptBadOption(optCon, 0), poptStrerror(c));
    return 1;
  }
  poptFreeContext(optCon);

  if (talk) {
    printf("Generating data for polyhedra talk figures.\n");
    return generate_talk_figs();
  }
  // ---------------------------------------------------------------------------
  // Verify we have reasonable arguments and set secondary parameters
  // ---------------------------------------------------------------------------
  const bool real_walls = !fake_walls;
  double len[3] = {periodic[0] + walls[0], periodic[1] + walls[1],
                         periodic[2] + walls[2]};

  const poly_shape shape(shape_name, ratio);
  if (shape.type == NONE) {
    fprintf(stderr, "\nInvalid shape.\n");
    return 1;
  }
  if(R == 0) {
    R = 1/(shape.vertices[1] - shape.vertices[0]).norm();
    printf("\nRadius was not specified, so setting radius to %g.\n", R);
  }
  if(ff != 0) {
    const double fac = R*pow(shape.volume*N/ff/len[0]/len[1]/len[2], 1.0/3.0);
    for(int i=0; i<3; i++) {
      periodic[i] *= fac;
      walls[i] *= fac;
      len[i] *= fac;
    }
    printf("\nFilling fraction was specified, so setting cell dimensions to (%g, %g, %g).\n", len[0], len[1], len[2]);

  }
  if (N <= 0 || iterations < 0 || R <= 0 || neighborR <= 0 || dr <= 0 || scale < 0 ||
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
  const double eta = (double)N*shape.volume*R*R*R/len[0]/len[1]/len[2];
  if (eta > 1) {
    fprintf(stderr, "\nYou're trying to cram too many polyhedra into the cell. They will never fit. Filling fraction: %g\n", eta);
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
  for(int i=1; i<argc; i++) {
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
  const int density_bins = round((len[0] + len[1] + len[2])/de_density);
  long *density_histogram = new long[density_bins]();

  const double min_len = min(len[0], min(len[1], len[2]));

  const int n_gfuncs = 4;
  const char *g_names[] = {"radial", "plebian", "special", "corner"};
  int n_cylinders[n_gfuncs] = {1, 0, 0, shape.nvertices};
  if(shape.type == CUBOID) {
    n_cylinders[1] = 4;
    n_cylinders[2] = 2;
  } else if(shape.type == TRUNCATED_TETRAHEDRON) {
    n_cylinders[1] = 4;
    n_cylinders[2] = 4;
  } else if(shape.type == TETRAHEDRON) {
    n_cylinders[1] = 4;
    n_cylinders[2] = 0;
  } else {
    n_cylinders[1] = shape.nfaces*2;
    n_cylinders[2] = 0;
  }

  const int g_bins = round(min_len/2/de_g);
  long *g_histogram = new long[n_gfuncs*g_bins]();

  // // this order parameter is the absolute value of the dot product of normal vector
  // // to the side most aligned with the z-axis with zhat.
  // const int order_bins = len[2]/de_density;
  // const double dcostheta = 1.0/order_bins;

  // long *order_parameter_histogram = new long[(int)round(len[2]/de_density)*order_bins]();


  polyhedron *polyhedra = new polyhedron[N];
  // long dZ = 0;

  // Initialize the random number generator with our seed
  random::seed(seed);

  // ---------------------------------------------------------------------------
  // Set up the initial grid of polyhedra
  // ---------------------------------------------------------------------------
  // Approximate the maximum number of neighbors a polyhedron could have
  // hokey guess, but it should always be high enough
  const double neighbor_sphere_vol = 4.0/3.0*M_PI*uipow(2.5*R+neighborR, 3) - shape.volume;
  int max_neighbors = min(N, neighbor_sphere_vol / shape.volume);

  for(int i=0; i<N; i++) {
    polyhedra[i].mypoly = &shape;
    polyhedra[i].R = R;
  }
  int max_attempts = 1;
  // cube
  if (shape.type == CUBE) {
    sprintf(structure, "cubic");
    int nx = ceil(pow((double)N*len[0]*len[0]/len[1]/len[2], 1.0/3.0));
    int ny = ceil((double)len[1]/len[0]*nx);
    int nz = ceil((double)len[2]/len[0]*nx);
    while (nx*ny*nz < N) {
      nx ++;
      ny ++;
      nz ++;
    }
    const double xspace = len[0]/double(nx);
    const double yspace = len[1]/double(ny);
    const double zspace = len[2]/double(nz);
    printf("Setting up grid with polyhedra: (%i, %i, %i) and space: (%g, %g, %g)\n",
           nx, ny, nz, xspace, yspace, zspace);

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
    }
    delete[] spot_used;
  }
  // non-cube cuboids
  else if (shape.nvertices == 8) {
    sprintf(structure, "cubic");
    const double xlen = fabs(2*shape.vertices[0][0]);
    const double ylen = fabs(2*shape.vertices[0][1]);
    const double zlen = fabs(2*shape.vertices[0][2]);

    const double ratio = zlen/len[2]/xlen*len[0];

    // want N_adj in the form a*a*a*ratio
    const int N_adjusted = round(ratio*ratio*uipow(ceil(pow(N/ratio/ratio, 1.0/3.0)), 3));

    int nx = round(pow(N_adjusted/ratio/ratio, 1.0/3.0)*ratio);
    int ny = round(pow(N_adjusted/ratio/ratio, 1.0/3.0)*ratio);
    int nz = round(pow(N_adjusted/ratio/ratio, 1.0/3.0));

    const double xspace = len[0]/nx;
    const double yspace = len[1]/ny;
    const double zspace = len[2]/nz;
    bool *spot_used = new bool[nx*ny*nz]();

    printf("%i, %i, %i, %i, %g, %g, %g\n", N_adjusted, nx, ny, nz, xspace, yspace, zspace);

    for(int i=0; i<N; i++) {
      int x, y, z;
      do {
        x = floor(random::ran()*nx);
        y = floor(random::ran()*ny);
        z = floor(random::ran()*nz);
      } while (spot_used[x*ny*nz + y*nz + z]);
      spot_used[x*ny*nz + y*nz + z] = true;
      polyhedra[i].pos = vector3d(x*xspace + xlen/2, y*yspace + ylen/2, z*zspace + zlen/2);
    }
  }
  // truncated tetrahedra
  else if (shape.type == TRUNCATED_TETRAHEDRON) {
    if (strcmp(structure, "arsenic") == 0) {
      // form a lattice that will be able to give us the maximum filling fraction
      // which is 207/208
      // Based on the paper Crystalline Assemblies... by Damasceno et al.

      const double fac = R/sqrt(11.0)/3.0*1.1;
      // lattice vectors:
      const vector3d e1 = vector3d(-4, 8, 12)*fac;
      const vector3d e2 = vector3d(12, -4, 8)*fac;
      const vector3d e3 = vector3d(8, 12, -4)*fac;

      const vector3d offset = vector3d(6, 6, 6)*fac;

      int nx = 2*ceil(pow((double)N*len[0]*len[0]/len[1]/len[2], 1.0/3.0));
      int ny = 2*ceil((double)len[1]/len[0]*nx);
      int nz = 2*ceil((double)len[2]/len[0]*nx);

      bool *spot_used = new bool[nx*ny*nz]();
      int x=0, y=0, z=0;
      for(int i=0; i<N; i+=2) {
        do {
          do {
            x = floor(random::ran()*nx);
            y = floor(random::ran()*ny);
            z = floor(random::ran()*nz);
          } while (spot_used[x*ny*nz + y*nz + z]);
          spot_used[x*ny*nz + y*nz + z] = true;

          x -= nx/2;
          y -= ny/2;
          z -= nz/2;
          polyhedra[i].pos = offset + x*e1 + y*e2 + z*e3;
          if (i < N-1) {
            polyhedra[i+1].pos = x*e1 + y*e2 + z*e3;
            polyhedra[i+1].rot = rotation(M_PI, vector3d(1, 1, 0));
          }
        } while((x*e1 + y*e2 + z*e3 + offset)[0] + R > len[0] || (x*e1 + y*e2 + z*e3 + offset)[1] + R > len[1] || (x*e1 + y*e2 + z*e3 + offset)[2] + R > len[2] || (x*e1 + y*e2 + z*e3)[0] < 0 || (x*e1 + y*e2 + z*e3)[1] < 0 || (x*e1 + y*e2 + z*e3)[2] - R < 0);
      }
      delete[] spot_used;
    }
    else if (strcmp(structure, "ice") == 0) {
      const double fac = R/sqrt(11.0);
      // lattice vectors:
      vector3d e1 = vector3d(0, 4, 4)*fac;
      vector3d e2 = vector3d(4, 0, 4)*fac;
      vector3d e3 = vector3d(4, 4, 0)*fac;


      // now we want to resize the lattice vectors so we'll have even spacing
      const double ne1_doub = sqrt(sqr(len[1]) + sqr(len[2]))/e1.norm();
      const int ne1 = int(ne1_doub);
      const double ne2_doub = sqrt(sqr(len[0]) + sqr(len[2]))/e2.norm();
      const int ne2 = int(ne2_doub);
      const double ne3_doub = sqrt(sqr(len[0]) + sqr(len[1]))/e3.norm();
      const int ne3 = int(ne3_doub);
      // N_adjusted is the smallest even cube that is at least as large as N
      const int N_adjusted = uipow(ceil(pow(N + N%2, 1.0/3.0)), 3);
      const double fraction = pow((double)N_adjusted/ne1/ne2/ne3, 1.0/3.0);

      const double space = (ne1_doub + ne2_doub + ne3_doub)/(ne1 + ne2 + ne3)/fraction;

      printf("%i, %i, %i\n", ne1, ne2, ne3);
      printf("%g, %g, %g\n", ne1_doub, ne2_doub, ne3_doub);
      e1 *= space;
      e2 *= space;
      e3 *= space;
      const vector3d offset = vector3d(2, 2, 2)*fac*space;
      bool *spot_used = new bool[8*ne1*ne2*ne3]();
      for(int i=0; i<N; i+=2) {
        int e1_i, e2_i, e3_i;
        vector3d pos;
        bool good_spot;
        int attempts = 0;
        do {
          e1_i = floor(random::ran()*2*ne1 - ne1);
          e2_i = floor(random::ran()*2*ne2 - ne2);
          e3_i = floor(random::ran()*2*ne3 - ne3);
          pos = e1_i*e1 + e2_i*e2 + e3_i*e3;
          good_spot = true;
          if (spot_used[(e1_i+ne1)*4*ne2*ne3 + (e2_i+ne2)*2*ne3 + e3_i+ne3])
            good_spot = false;
          else {
            for(int j=0; j<3; j++){
              if(pos[j] < 0 || pos[j] + offset[j] > len[j]) {
                good_spot = false;
                break;
              }
            }
          }
          attempts ++;
        } while(!good_spot);
        max_attempts = max(attempts, max_attempts);
        polyhedra[i].pos = pos;
        polyhedra[i].rot = rotation(M_PI, vector3d(1, 1, 0));
        if(i < N-1) polyhedra[i+1].pos = offset + pos;
        spot_used[(e1_i+ne1)*4*ne2*ne3 + (e2_i+ne2)*2*ne3 + e3_i+ne3] = true;
      }
      delete[] spot_used;

      for(int i=0; i<N; i++) {
        polyhedra[i].pos = fix_periodic(polyhedra[i].pos, len);
      }
    }
  }
  took("Placement");
  printf("The most attempts for a polyhedron was %i.\n", max_attempts);

  // ---------------------------------------------------------------------------
  // Save the initial configuration for troubleshooting
  // ---------------------------------------------------------------------------
  char *vertices_fname = new char[1024];
  sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape.name, N, -1);
  save_locations(polyhedra, N, vertices_fname, len);
  delete[] vertices_fname;

  int most_neighbors =
    initialize_neighbor_tables(polyhedra, N, neighborR + 2*dr, max_neighbors, periodic);
  if (most_neighbors < 0) {
    fprintf(stderr, "The guess of %i max neighbors was too low. Exiting.\n", max_neighbors);
    return 1;
  }
  printf("Neighbor tables initialized. The most neighbors is %i, whereas the max allowed is %i.\n", most_neighbors, max_neighbors);

  print_all(polyhedra, N, periodic);
  print_bad(polyhedra, N, periodic);

  // ---------------------------------------------------------------------------
  // Make sure no polyhedra are overlapping
  // ---------------------------------------------------------------------------
  for(int i=0; i<N; i++) {
    if (!in_cell(polyhedra[i], walls, real_walls)) {
      printf("Oops, this is embarassing. I seem to have placed some things outside our cell.\n");
      printf("You might want to look into that.\n");
      return 17;
    }
    for(int j=i+1; j<N; j++) {
      if (overlap(polyhedra[i], polyhedra[j], periodic)) {
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
    // Move each polyhedron once
    // ---------------------------------------------------------------
    for(int i=0; i<N; i++) {
      totalmoves ++;
      int move_val = move_one_polyhedron(i, polyhedra, N, periodic, walls, real_walls,
                                        neighborR, scale, theta_scale, max_neighbors, dr);
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
      printf("Iteration %li, acceptance rate of %g, scale: %g.\n", iteration, (double)workingmoves/totalmoves, scale);
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
      fflush(stdout);
    }
  }
  took("Initialization");

  // ---------------------------------------------------------------------------
  // Save the post-initialization configuration for troubleshooting
  // ---------------------------------------------------------------------------
  vertices_fname = new char[1024];
  sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape.name, N, -2);
  save_locations(polyhedra, N, vertices_fname, len);
  delete[] vertices_fname;

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
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*30; // top out at half hour interval
  clock_t last_output = clock(); // when we last output data

  int frame = 0;
  totalmoves = 0, workingmoves = 0, old_totalmoves = 0, old_workingmoves = 0;
  neighbor_updates = 0, neighbor_informs = 0;

  for(long iteration=1; iteration<=iterations; iteration++) {
    // ---------------------------------------------------------------
    // Move each polyhedron once
    // ---------------------------------------------------------------
    for(int i=0; i<N; i++) {
      totalmoves ++;
      int move_val = move_one_polyhedron(i, polyhedra, N, periodic, walls, real_walls,
                                         neighborR, scale, theta_scale, max_neighbors, dr);
      workingmoves += move_val & 1;
    }
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

      // // Order parameter:
      // if (shape.type == CUBE) {
      //   double costheta = 0;
      //   for(int j=0; j<polyhedra[i].mypoly->nfaces; j++) {
      //     const double new_costheta = polyhedra[i].rot.rotate_vector(polyhedra[i].mypoly->faces[j])[2];
      //     costheta = max(costheta, fabs(new_costheta));
      //   }
      //   const int costheta_i = costheta/dcostheta;
      //   order_parameter_histogram[z_i*order_bins + costheta_i] ++;
      // }
    }
    // Distribution histogram:
    // g_names = {"radial", "plebian", "special", "corner"}

    // unrot lets us work with treating the first polyhedron as sitting at the origin,
    // unrotated

    // this is an O(N^2) calculation, so we're only going to do it every N^2 moves
    if(iteration %N == 0) {
      for(int i=0; i<N; i++) {
        const rotation unrot = polyhedra[i].rot.conj();
        for(int j=0; j<N; j++) {
          if(i != j) {
            // note: be careful about orientation when adding that
            const vector3d pos2 = unrot.rotate_vector
              (periodic_diff(polyhedra[i].pos, polyhedra[j].pos, periodic));
            // radial
            const double r = pos2.norm();
            const int r_i = floor(r/de_g);
            if(r_i < g_bins) g_histogram[r_i] ++;

            // plebian face
            // for cuboids, this is +-xhat and +-yhat
            // for truncated_tetrahedra, this is the face vectors
            // for everythig else, this is all the faces
            if(shape.type == CUBOID) {
              const vector3d faces[] = {vector3d(1,0,0), vector3d(0,1,0)};
              for(int k=0; k<2; k++) {
                const double e = pos2[k];
                const double dist = (pos2 - e*faces[k]).norm();
                if (dist < dr_g) {
                  const int e_i = floor(fabs(e/de_g));
                  if (e_i < g_bins) g_histogram[g_bins + e_i] ++;
                }
              }
            } else {
              for(int k=0; k<shape.nfaces; k++) {
                double e = pos2.dot(shape.faces[k]);
                const double dist = (pos2 - e*shape.faces[k]).norm();
                if (dist < dr_g) {
                  if(shape.type != TETRAHEDRON && shape.type != TRUNCATED_TETRAHEDRON)
                    e = fabs(e);
                  const int e_i = floor(e/de_g);
                  if(e_i >= 0 && e_i < g_bins) g_histogram[g_bins + e_i] ++;
                }
              }
            }
            // special face
            // for cuboids, this is +- zhat
            // for truncated tetrahedra, this is the negative of the face vectors
            // for everything else, this is nothing
            if(shape.type == CUBOID) {
              const double z = pos2[2];
              const double zdist = (pos2 - z*vector3d(0,0,1)).norm();
              if (zdist < dr_g) {
                const int z_i = floor(fabs(z/de_g));
                if(z_i < g_bins) g_histogram[2*g_bins + z_i] ++;
              }
            } else if (shape.type == TRUNCATED_TETRAHEDRON) {
              for(int k=0; k<shape.nfaces; k++) {
                double e = pos2.dot(-shape.faces[k]);
                const double dist = (pos2 + e*shape.faces[k]).norm();
                if (dist < dr_g) {
                  const int e_i = floor(e/de_g);
                  if(e_i >= 0 && e_i < g_bins) g_histogram[2*g_bins + e_i] ++;
                }
              }
            }
            // vertices
            for(int k=0; k<shape.nvertices; k++) {
              const double e = pos2.dot(shape.vertices[k]);
              const double dist = (pos2 - e*shape.vertices[k]).norm();
              if (dist < dr_g) {
                const int e_i = floor(e/de_g);
                if(e_i >= 0 && e_i < g_bins) g_histogram[3*g_bins + e_i] ++;
              }
            }
          }
        }
      }
    }
    // ---------------------------------------------------------------
    // Get pressure info - this is slow, might want to do it less often
    // fixme - disabled for now
    // ---------------------------------------------------------------
    // for(int i=0; i<N; i++) {
    //   dZ += overlaps_with_any(polyhedra[i], polyhedra, periodic, true, dr);
    // }
    // dZ /= 2; // since we check each thing against each of its neighbors, we end up
    //          // double counting. Ideally, we would count half as much instead of dividing
    // ---------------------------------------------------------------
    // Save to file
    // ---------------------------------------------------------------
    const clock_t now = clock();
    if ((now - last_output > output_period) || iteration==iterations) {
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
      sprintf(density_fname, "%s/%s-density-%s-%i.dat", dir, filename, shape.name, N);
      FILE *densityout = fopen((const char *)density_fname, "w");
      delete[] density_fname;
      fprintf(densityout, "%s", headerinfo);
      fprintf(densityout, "%s", countinfo);
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

      // Save distribution funtions
      // fixme: assumes homogeneous for density
      char *g_fname = new char[1024];
      sprintf(g_fname, "%s/%s-g-%s-%i.dat", dir, filename, shape.name, N);
      FILE *g_out = fopen((const char *)g_fname, "w");
      delete[] g_fname;
      fprintf(g_out, "%s", headerinfo);
      fprintf(g_out, "%s", countinfo);
      fprintf(g_out, "\ne       ");
      const double density = N/len[0]/len[1]/len[2];
      const double vol0 = len[0]*len[1]*len[2];
      for(int i=0; i<n_gfuncs; i++) fprintf(g_out, "%8s  ", g_names[i]);
      fprintf(g_out, "\n");
      for(int e_i=0; e_i<g_bins; e_i++) {
        const double e = (e_i + 0.5)*de_g;
        fprintf(g_out, "%6.3f  ", e);
        for(int i=0; i<n_gfuncs; i++) {
          const double vol1 = i==0 ? 4.0/3.0*M_PI*(uipow(e+de_g/2, 3) - uipow(e-de_g/2, 3))
                                   : M_PI*sqr(dr_g)*de_g;
          const double probability = (double)g_histogram[i*g_bins + e_i]/totalmoves;
          const double n2 = probability/vol0/vol1;
          const double g = n2/density/density*N*N/n_cylinders[i];
          fprintf(g_out, "%8.5f  ", g);
        }
        fprintf(g_out, "\n");
      }
      fclose(g_out);
      // // Save order paramter
      // if (shape.type == CUBE) {
      //   char *order_fname = new char[1024];
      //   sprintf(order_fname, "%s/%s-order-%s-%i.dat", dir, filename, shape.name, N);
      //   FILE *order_out = fopen((const char *)order_fname, "w");
      //   delete[] order_fname;
      //   fprintf(order_out, "%s", headerinfo);
      //   fprintf(order_out, "%s", countinfo);
      //   for(int z_i=0; z_i<round(len[2]/de_density); z_i++) {
      //     for(int costheta_i=0; costheta_i<order_bins; costheta_i++) {
      //       const double costheta = (double)order_parameter_histogram[z_i*order_bins + costheta_i]*N/density_histogram[xbins + ybins + z_i]/dcostheta;
      //       fprintf(order_out, "%4.2f ", costheta);
      //     }
      //     fprintf(order_out, "\n");
      //   }
      //   fclose(order_out);
      // }

      // // Save pressure
      // char *pressure_fname = new char[1024];
      // sprintf(pressure_fname, "%s/%s-pressure-%s-%i.dat", dir, filename, shape.name, N);
      // FILE *pressureout = fopen((const char *)pressure_fname, "a");
      // delete[] pressure_fname;
      // const double dV = N*shape.volume*(uipow(R+dr, 3) - uipow(R, 3));
      // const double pressure = dZ/dV/totalmoves;
      // fprintf(pressureout, "%li %g %li\n", totalmoves, pressure, dZ);
      // old_totalmoves = totalmoves;
      // dZ = 0;
      // fclose(pressureout);

      delete[] countinfo;
      // ---------------------------------------------------------------------------
      // Save the configuration for troubleshooting
      // ---------------------------------------------------------------------------
      vertices_fname = new char[1024];
      sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape.name, N, -3);
      save_locations(polyhedra, N, vertices_fname, len);
      delete[] vertices_fname;
    }
    // ---------------------------------------------------------------
    // Save locations of vertices if desired
    // ---------------------------------------------------------------
    if (vertex_period > 0 && iteration % vertex_period == 0) {
      printf("Saving vertex locations. Frame: %i. Iteration: %li.\n", frame, iteration);
      vertices_fname = new char[1024];
      sprintf(vertices_fname, "%s/vertices/%s-vertices-%s-%i-%i.dat", dir, filename, shape.name, N, frame);
      char *comment = new char[1024];
      sprintf(comment, "period: %i, iteration: %li", vertex_period, iteration);
      save_locations(polyhedra, N, vertices_fname, len, comment);
      delete[] vertices_fname;
      frame ++;
    }
    // ---------------------------------------------------------------
    // Print out timing information if desired
    // ---------------------------------------------------------------
    if (totime > 0 && iteration % totime == 0) {
      char *iter = new char[1024];
      sprintf(iter, "%i iterations", totime);
      took(iter);
      delete[] iter;
      printf("Iteration %li, acceptance rate of %g, scale: %g.\n", iteration, (double)workingmoves/totalmoves, scale);
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
      fflush(stdout);
    }
  }
  // ---------------------------------------------------------------------------
  // END OF MAIN PROGRAM LOOP
  // ---------------------------------------------------------------------------
  print_bad(polyhedra, N, periodic);

  //delete[] polyhedra; fixme
  delete[] density_histogram;
  delete[] g_histogram;

  delete[] headerinfo;
  return 0;
}
// -----------------------------------------------------------------------------
// END OF MAIN
// -----------------------------------------------------------------------------

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

void save_locations(const polyhedron *p, int N, const char *fname, const double len[3], const char *comment) {
  FILE *out = fopen((const char *)fname, "w");
  fprintf(out, "# %s\n", comment);
  fprintf(out, "%g %g %g\n", len[0], len[1], len[2]);
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

int generate_talk_figs() {
  const poly_shape cube("cube");
  const poly_shape truncated_tetrahedron("truncated_tetrahedron");
  const poly_shape tetrahedron("tetrahedron");

  // Generate background image for title slide
  int N = 7;
  polyhedron *p = new polyhedron[N];
  const double len[3] = {20, 20, 20};
  const double R = 1;
  for(int i=0; i<N; i++) {
    p[i].R = R;
  }
  p[0].mypoly = &truncated_tetrahedron;
  p[1].mypoly = &cube;
  p[2].mypoly = &truncated_tetrahedron;
  p[3].mypoly = &tetrahedron;
  p[4].mypoly = &cube;
  p[5].mypoly = &truncated_tetrahedron;
  p[6].mypoly = &tetrahedron;

  const double periodic[3] = {0, 0, 0};
  const double walls[3] = {6.5, 4, 5};
  const bool real_walls = true;
  const int max_neighbors = N;
  const double scale = 1.5;
  const double theta_scale = M_PI/4.0;
  const double neighborR = .5;

  for(int i=0; i<N; i++) {
    p[i].pos = vector3d(random::ran(), random::ran(), random::ran())*3;
  }
  initialize_neighbor_tables(p, N, neighborR, max_neighbors, periodic);

  for(int j=0; j<100; j++) {
    for(int i=0; i<N; i++) {
      move_one_polyhedron(i, p, N, periodic, walls, real_walls,
                          neighborR, scale, theta_scale, max_neighbors, 0);
    }
  }
  for(int i=0; i<N; i++) {
    p[i].pos += vector3d(0, -4.5, 0);
  }
  save_locations(p, N, "talks/polyhedra/dat/background.dat", len);

  delete[] p;
  // Generate image of ice structure

  N = 26;
  p = new polyhedron[N];


  for(int i=0; i<N; i++) {
    p[i].R = R;
    p[i].mypoly = &truncated_tetrahedron;
  }
  const double fac = R/sqrt(11.0)*1.05;
  // lattice vectors:
  const vector3d e1 = vector3d(0, 4, 4)*fac;
  const vector3d e2 = vector3d(4, 0, 4)*fac;
  const vector3d e3 = vector3d(4, 4, 0)*fac;
  const vector3d offset = vector3d(2, 2, 2)*fac;

  double rad = 3.5;
  int nx=10, ny=10, nz=10;
  int x=-nx/2, y=-ny/2, z=-nz/2;
  int i=0;
  while(i < N-1) {
    if((x*e1 + y*e2 + z*e3).norm() < rad) {
      p[i].pos = x*e1 + y*e2 + z*e3;
      p[i].rot = rotation(M_PI, vector3d(1, 1, 0));
      if(i < N-1) p[i+1].pos = offset + x*e1 + y*e2 + z*e3;
      else break;
      i += 2;
    }
    x++;
    if(x >= nx/2) {
      x = -nx/2;
      y++;
      if(y >= ny/2) {
        y = -ny/2;
        z++;
        if(z >= nz/2) break;
      }
    }
  }
  save_locations(p, N, "talks/polyhedra/dat/ice-structure.dat", len);

  delete[] p;

  // tetrahedra images
  p = new polyhedron[1];
  rotation *rots = new rotation[3];
  rots[1] = rotation(M_PI/4, vector3d(0, 1, 0));
  rots[2] = rotation(M_PI/3, vector3d(1, 0, 0));
  for(int i=0; i<3; i++) {
    p[0].R = R;
    p[0].mypoly = &truncated_tetrahedron;
    p[0].rot = rots[i];
    char *fname = new char[1024];
    sprintf(fname, "talks/polyhedra/dat/tet-%i.dat", i);
    save_locations(p, 1, fname, len);
  }

  return 0;
}
