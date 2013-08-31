#include <stdio.h>
#include <time.h>
#include "MersenneTwister.h"
#include "vector3d.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utilities.h"

struct poly_shape {
  int nvertices;
  int nfaces;
  vector3d vertices[20];
  vector3d faces[6];
  char name[1024];
};
// Note: faces is a list of normal vectors perpendicular to the actual faces.
// The sign of this vector is unimportant, so any two parallel faces are only
// counted as one face together.

struct polyhedron {
  vector3d pos;
  quaternion rot;
  double R;
  poly_shape *mypoly;
};

// Global Constants
const bool debug = false; // prints out a lot of extra information
                         // also performs extra computation
                         // should be false except when something is wrong

// Global "Constants" -- set at runtime then unchanged
long unsigned int seed = 0; // for random number generator
int iterations;
int N;
bool periodic[3] = {false, false, false};
bool walls[3] = {false, false, false};
bool real_walls = true; // true means shapes can't overlap wall at all,
                        // false means just centers can't overlap wall
double len[3] = {20, 20, 20};
double dw_density = 0.1;
double R = 1.0;
double scale = 0.005;
double theta_scale = 0.005;


// All the shapes:
poly_shape cube;
poly_shape tetrahedron;
poly_shape truncated_tetrahedron;

// Functions:

inline double min(double a, double b) { return (a < b)? a : b; }
inline double max(double a, double b) { return (a > b)? a : b; }

// Sets up the vertices, faces, etc. for each polyhedron. Only run once.
void define_shapes();

// Check whether two polyhedra overlap
inline bool overlap(const polyhedron &a, const polyhedron &b);

// Check whether a polyhedron overlaps with any in the array, except for
// the nth one. Useful for post-move testing with a temporary polyhedron.
inline bool overlaps_with_any(const polyhedron &a, const polyhedron *bs, int n);

// Check if polyhedron is outside its allowed locations
inline bool in_cell(const polyhedron &p);

// Return the vector pointing from a to b, accounting for periodic
// boundaries
inline vector3d periodic_diff(const vector3d &a, const vector3d  &b);

// Move v in a random direction by a distance determined by a gaussian distribution
inline vector3d move(const vector3d &v, double scale);

// Rotate q in a random direction by an amount determined by a gaussian distribution
inline quaternion rotate(const quaternion &q, double scale);

// Rotate the vector v about the axis of vector part of q
// by an amount equal to the scalar part of q
inline vector3d rotate_vector(const vector3d &v, const quaternion &q);

// Generate a random number in the range [0, 1) using a fixed seed
inline double ran();

// Generate a random point in a gaussian distribution
inline vector3d ran3();

// If v is outside the cell, and there are periodic boundary condition(s), it is
// moved (appropriately) into the cell
inline vector3d fix_periodic(vector3d newv);

// States how long it's been since last took call.
static void took(const char *name);



// The following functions only do anything if debug is true:

// Prints the locations, rotations, radii, and shapes of every polyhedron
// As well as if any overlap or are outside the cell.
inline void print_all(const polyhedron *p);

// Only print those shapes that overlap or are outside the cell
// Also prints those they overlap with
inline void print_bad(const polyhedron *p);




int main(int argc, const char *argv[]) {
  define_shapes();

  poly_shape *shape;
  shape = &cube;
  // input parameters
  const int minparams = 4;
  if (argc < minparams) {
    printf("\nUse: %s N iterations filename\n", argv[0]);
    printf("Additional optional parameters:\n\t\tperiodx, periody, periodz, wallx, wally, wallz,\n");
    printf("\t\tdimensions {LENX LENY LENZ}, dw_density {VAL}, seed {VAL}, shape {POLY},\n");
    printf("\t\tR {VAL}, scale {VAL}, theta_scale {VAL}, fake_walls\n");
    printf("Anything in brackets means to put values there.\n");
    printf("Available shapes currently include: cube.\n");
    return 1;
  }
  N = atoi(argv[1]);
  const long iterations = atol(argv[2]);
  const char *filename = argv[3];
  printf("----------------------------------------------------------------------\n");
  printf("Running %s with the following parameters:\n", argv[0]);
  printf("Setting number of polyhedra to %i\n", N);
  printf("Running for %4.2e iterations\n", double(iterations));
  printf("Saving data to: %s\n", filename);
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
    } else if (strcmp(argv[i], "dw_density") == 0) {
      dw_density = atof(argv[i+1]);
      printf("Setting density resolution to %g.\n", dw_density);
      i += 1;
    } else if (strcmp(argv[i], "R") == 0) {
      R = atof(argv[i+1]);
      printf("Setting radius of the circumscribed sphere to %g.\n", R);
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
    } else if (strcmp(argv[i], "shape") == 0) {
      if (strcmp(argv[i+1], "cube") == 0) {
        shape = &cube;
        printf("Using cubes.\n");
      } else if (strcmp(argv[i+1], "tetrahedron") == 0) {
        shape = &tetrahedron;
        printf("Using cubes.\n");
      } else {
        printf("Invalid shape: %s\n.", argv[i+1]);
        return 1;
      }
      i += 1;
    } else {
      printf("Invalid parameter: %s.\n", argv[i]);
      return 1;
    }
  }
  if (debug) printf("DEBUG MODE IS ENABLED!\n");
  else printf("Debug mode disabled.\n");
  printf("----------------------------------------------------------------------\n\n");
  const int density_bins = len[2]/dw_density;
  long *density_histogram = new long[density_bins]();
  polyhedron *polyhedra = new polyhedron[N];
  polyhedron initialpolys[N];
  double displacement[N];
  long totalmoves = 0, workingmoves = 0;

  ////////////////////////////////////////////////////////////////////////////////
  // Start by arranging the shapes in a grid.
  // Right now is set up to simply make them evenly spaced with no rotation.
  // It will fail with an error if any of them overlap.
  // NOTE: This only really works if the cell is at least nearly square,
  // and if there's enough room to fit the shapes without rotation.
  // A better version, that is shape-dependent, will be required for higher
  // filling fractions.
  const int nx = ceil(pow(N*len[0]*len[0]/len[1]/len[2], 1.0/3.0));
  const int ny = ceil((len[1]/len[0])*nx);
  const int nz = ceil(N/double(nx*ny));
  const double xspace = len[0]/double(nx);
  const double yspace = len[1]/double(ny);
  const double zspace = len[2]/double(nz);
  if(debug) printf("n polyhdra: (%i, %i, %i), with space: (%g, %g, %g)\n",
                   nx, ny, nz, xspace, yspace, zspace);
  int x = 0, y = 0, z = 0;
  for(int i=0; i<N; i++) {
    polyhedra[i].mypoly = shape;
    polyhedra[i].R = R;
    polyhedra[i].pos[0] = (x + 0.5)*xspace;
    polyhedra[i].pos[1] = (y + 0.5)*yspace;
    polyhedra[i].pos[2] = (z + 0.5)*zspace;

    initialpolys[i].mypoly = shape;
    initialpolys[i].R = R;
    initialpolys[i].pos[0] = (x + 0.5)*xspace;
    initialpolys[i].pos[1] = (y + 0.5)*yspace;
    initialpolys[i].pos[2] = (z + 0.5)*zspace;

    x ++;
    if (x >= nx) {
      x = 0; y ++;
      if (y >= ny) {
        y = 0; z ++;
      }
    }
  }
  print_all(polyhedra);
  print_bad(polyhedra);
  took("Initial setup");
  // Make sure none are overlapping:
  for(int i=0; i<N; i++) {
    if (!in_cell(polyhedra[i])) {
      printf("Oops, this is embarassing. I seem to have placed some things outside our cell.\n");
      printf("You might want to look into that.\n");
      return 1;
    }
    for(int j=i+1; j<N; j++) {
      if (overlap(polyhedra[i], polyhedra[j])) {
        printf("ERROR in intial placement. We have overlaps!!!\n");
        printf("AHHHHHH I DON'T KNOW WHAT TO DO!@!!!!1111\n");
        return 1;
      }
    }
  }
  fflush(stdout);
  // fixme: multiply output_period by 60
  clock_t output_period = CLOCKS_PER_SEC; // start at outputting every minute
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*60; // top out at one hour interval
  clock_t last_output = clock(); // when we last output data
  bool initial_phase = true;
  long start_keeping = 0;
  // BEGIN MAIN PROGRAM LOOP ////////////////////////////////////////////////////////
  for(long iteration=1; iteration<=iterations; iteration++) {
    if (debug) {
      print_bad(polyhedra);
      if (iteration%10000 == 0) {
        print_all(polyhedra);
        double maxdis = 0, mindis = 1000000, sumdis = 0;
        for(int i=0; i<N; i++) {
          maxdis = max(maxdis, displacement[i]);
          mindis = min(mindis, displacement[i]);
          sumdis += displacement[i];
        }
        printf("biggest displacement: %4.2f, smallest: %4.2f, mean: %4.2f\n", maxdis, mindis, sumdis/N);
      }
      if (iteration%1000 == 0) {
        char *iter = new char[1024];
        sprintf(iter, "One thousand iterations (current iteration: %li)", iteration);
        took(iter);
        delete[] iter;
        fflush(stdout);
      }
    }

    // MOVING SHAPES ------------------------------------------------------
    for(int i=0; i<N; i++) {
      totalmoves++;
      polyhedron temp;
      temp.pos = move(polyhedra[i].pos, scale);
      temp.rot = rotate(polyhedra[i].rot, theta_scale);
      temp.R = polyhedra[i].R;
      temp.mypoly = polyhedra[i].mypoly;
      if (in_cell(temp)) {
        const bool overlaps = overlaps_with_any(temp, polyhedra, i);
        if (!overlaps) {
          polyhedra[i] = temp;
          workingmoves ++;
          if (debug)
            displacement[i] = (polyhedra[i].pos-initialpolys[i].pos).norm();
        }
      }
    }
    // STOP and go to next iteration if we don't want to store data yet
    // Since more start at lower z-values, once there are more at higher values
    // We should be ready
    if (initial_phase) {
      if (iteration%1000 == 0) {
        int count1 = 0, count2 = 0;
        for(int i=0; i<N; i++) {
          if (polyhedra[i].pos[2] < len[2]/2.0)
            count1 ++;
          else
            count2 ++;
        }
        if (count1 >= count2) {
          printf("Not ready to start storing data yet - c1: %i, c2: %i, iteration: %li, acceptance rate: %4.2f\n", count1, count2, iteration, double(workingmoves)/totalmoves);
        }
        else {
          initial_phase = false;
          start_keeping = iteration;
          printf("\nTime to start data collection! c1: %i, c2: %i, iteration: %li\n", count1, count2, iteration);
          took("Initial movements");
          printf("\n");
        }
        fflush(stdout);
      }
    }
    else {
      // ADDING DATA TO HISTOGRAM(S) ----------------------------------------
      for(int i=0; i<N; i++) {
        // Density histogram:
        const int z_i = floor(polyhedra[i].pos[2]/dw_density);
        density_histogram[z_i] ++;
      }
      // SAVING TO FILE -----------------------------------------------------
      const clock_t now = clock();
      if ((now > last_output + output_period) || iteration==iterations) {
        last_output = now;
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

        const long count = N*(iteration - start_keeping);
        char *headerinfo = new char[4096];
        sprintf(headerinfo, "# dim: (%5.2f, %5.2f, %5.2f), period: (%i, %i, %i), \
walls: (%i, %i, %i), dw: %g, seed: %li\n# R: %f, scale: %g, theta_scale: %g, \
real_walls: %i\n# num counts: %li, iteration: %li, \
workingmoves: %li, totalmoves: %li, acceptance rate: %g\n",
                len[0], len[1], len[2], periodic[0], periodic[1], periodic[2],
                walls[0], walls[1], walls[2], dw_density, seed, R,
                scale, theta_scale, real_walls, count,
                iteration, workingmoves, totalmoves, double(workingmoves)/totalmoves);
        // saving density
        char *density_filename = new char[1024];
        sprintf(density_filename, "%s-density-%s-%i.dat", filename, shape->name, N);
        FILE *densityout = fopen((const char *)density_filename, "w");
        fprintf(densityout, "%s", headerinfo);
        fprintf(densityout, "\n# z     density   histogram\n");
        for(int z_i = 0; z_i < density_bins; z_i ++) {
          const double z = (z_i + 0.5)*dw_density;
          const double shell_volume = len[0]*len[1]*dw_density;
          const double density = (double)density_histogram[z_i]*N/count/shell_volume;
          fprintf(densityout, "%5.2f   %07.5f   %li\n", z, density, density_histogram[z_i]);
        }
        fclose(densityout);

        delete[] headerinfo;
      }
    }
  }
  print_bad(polyhedra);
  // END MAIN PROGRAM LOOP //////////////////////////////////////////////////////////

  delete[] polyhedra;
  delete[] density_histogram;
  return 0;
}
// END OF MAIN **********************************************************************







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
    const vector3d axis = rotate_vector(a.mypoly->faces[i], a.rot);
    double projection = axis.dot(rotate_vector(a.mypoly->vertices[0]*a.R, a.rot));
    double mina = projection, maxa = projection;
    for (int j=1; j<a.mypoly->nvertices; j++) {
      projection = axis.dot(rotate_vector(a.mypoly->vertices[j]*a.R, a.rot));
      if (projection < mina) mina = projection;
      else if (projection > maxa) maxa = projection;
    }
    projection = axis.dot(rotate_vector(b.mypoly->vertices[0]*b.R, b.rot) + ab);
    double minb = projection, maxb = projection;
    for (int j=1; j<b.mypoly->nvertices; j++) {
      projection = axis.dot(rotate_vector(b.mypoly->vertices[j]*b.R, b.rot) + ab);
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
    const vector3d axis = rotate_vector(b.mypoly->faces[i], b.rot);
    double projection = axis.dot(rotate_vector(a.mypoly->vertices[0]*a.R, a.rot));
    double mina = projection, maxa = projection;
    for (int j=1; j<a.mypoly->nvertices; j++) {
      projection = axis.dot(rotate_vector(a.mypoly->vertices[j]*a.R, a.rot));
      if (projection < mina) mina = projection;
      else if (projection > maxa) maxa = projection;
    }
    projection = axis.dot(rotate_vector(b.mypoly->vertices[0]*b.R, b.rot) + ab);
    double minb = projection, maxb = projection;
    for (int j=1; j<b.mypoly->nvertices; j++) {
      projection = axis.dot(rotate_vector(b.mypoly->vertices[j]*b.R, b.rot) + ab);
      if (projection < minb) minb = projection;
      else if (projection > maxb) maxb = projection;
    }
    if (mina > maxb || minb > maxa) {
      return false;
    }
  }
  return true;
}

inline bool overlaps_with_any(const polyhedron &a, const polyhedron *bs, int n) {
  // construct axes from a and a's projection onto them
  vector3d aaxes[a.mypoly->nfaces];
  double amins[a.mypoly->nfaces], amaxes[a.mypoly->nfaces];
  for (int i=0; i<a.mypoly->nfaces; i++) {
    aaxes[i] = rotate_vector(a.mypoly->faces[i], a.rot);
    double projection = aaxes[i].dot(rotate_vector(a.mypoly->vertices[0]*a.R, a.rot));
    amins[i] = projection, amaxes[i] = projection;
    for (int j=1; j< a.mypoly->nvertices; j++) {
      projection = aaxes[i].dot(rotate_vector(a.mypoly->vertices[j]*a.R, a.rot));
      amins[i] = min(projection, amins[i]);
      amaxes[i] = max(projection, amaxes[i]);
    }
  }
  for (int k=0; k<N; k++) {
    if (k != n) {
      const vector3d ab = periodic_diff(a.pos, bs[k].pos);
      if (ab.normsquared() < (a.R + bs[k].R)*(a.R + bs[k].R)) {
        bool overlap = true; // assume overlap until we prove otherwise or fail to
        // check projection of b against a's axes
        for (int i=0; i<a.mypoly->nfaces; i++) {
          double projection = aaxes[i].dot
            (rotate_vector(bs[k].mypoly->vertices[0]*bs[k].R, bs[k].rot) + ab);
          double bmin = projection, bmax = projection;
          for (int j=1; j<bs[k].mypoly->nvertices; j++) {
            projection = aaxes[i].dot
              (rotate_vector(bs[k].mypoly->vertices[j]*bs[k].R, bs[k].rot) + ab);
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
            const vector3d axis = rotate_vector(bs[k].mypoly->faces[i], bs[k].rot);
            double projection = axis.dot(rotate_vector(a.mypoly->vertices[0]*a.R, a.rot));
            double amin = projection, amax = projection;
            for (int j=1; j<a.mypoly->nvertices; j++) {
              projection = axis.dot(rotate_vector(a.mypoly->vertices[j]*a.R, a.rot));
              amin = min(projection, amin);
              amax = max(projection, amax);
            }
            projection = axis.dot
              (rotate_vector(bs[k].mypoly->vertices[0]*bs[k].R, bs[k].rot) + ab);
            double bmin = projection, bmax = projection;
            for (int j=1; j<bs[k].mypoly->nvertices; j++) {
              projection = axis.dot
                (rotate_vector(bs[k].mypoly->vertices[j]*bs[k].R, bs[k].rot) + ab);
              bmin = min(projection, bmin);
              bmax = max(projection, bmax);

            }
            if (amin > bmax || bmin > amax) {
              overlap = false;
              i = bs[k].mypoly->nfaces; //no overlap, move on to next
            }
          }
        }
        if (overlap) return true;
      }
    }
  }
  return false;
}

inline bool in_cell(const polyhedron &p) {
  if(real_walls) {
    for (int i=0; i<3; i++) {
      if (walls[i]) {
        if (p.pos[i] - p.R > 0.0 && p.pos[i] + p.R < len[i]) {
          continue;
        }
        double coord = (rotate_vector(p.mypoly->vertices[0]*p.R, p.rot) + p.pos)[i];
        double pmin = coord, pmax = coord;
        for (int j=1; j<p.mypoly->nvertices; j++) {
          coord = (rotate_vector(p.mypoly->vertices[j]*p.R, p.rot) + p.pos)[i];
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
  const vector3d newv = v + scale*ran3();
  return fix_periodic(newv);
}

inline quaternion rotate(const quaternion &q, double scale) {
  double x, y, z, r2, sintheta;
  do {
    do {
      x = 2*ran() - 1;
      y = 2*ran() - 1;
      r2 = x*x + y*y;
    } while (r2 >= 1 || r2 == 0);
    const double fac = sqrt(-2*log(r2)/r2);
    sintheta = fac*x*scale;
  } while (sintheta <= -1 or sintheta >= 1);
  const double costheta = sqrt(1 - sintheta*sintheta);
  do {
    x = 2*ran() - 1;
    y = 2*ran() - 1;
    z = 2*ran() - 1;
    r2 = x*x + y*y + z*z;
  } while(r2 >= 1 || r2 == 0);
  const double vfac = sintheta/sqrt(r2);
  const quaternion rot(costheta, x*vfac, y*vfac, z*vfac);
  if (debug) {
    quaternion totalrot = rot*q;
    if (fabs(1.0-rot.normsquared()) > 1e-15 || totalrot[0] > 1 || totalrot[0] < -1) {
      printf("%g\n", 1.0 - rot.normsquared());
      char rotstr[1024];
      rot.tostr(rotstr);
      double theta = 2*acos(rot[0]);
      printf("Rot: %s, theta: %4.2f, n: %.1f ", rotstr, theta, rot.norm());
      totalrot.tostr(rotstr);
      theta = 2*acos(totalrot[0]);
      printf("Total: %s, theta: %4.2f, n: %.1f\n", rotstr, theta, totalrot.norm());
      fflush(stdout);
    }
  }
  return (rot*q).normalized();
}

inline vector3d rotate_vector(const vector3d &v, const quaternion &q) {
  const quaternion product = q*quaternion(0, v)*q.conj();
  return vector3d(product[1], product[2], product[3]);
}

inline double ran() {
  static MTRand my_mtrand(seed); // always use the same random number generator (for debugging)!
  return my_mtrand.randExc(); // which is the range of [0,1)
}

inline vector3d ran3() {
  double x, y, r2;
  do {
    x = 2*ran() - 1;
    y = 2*ran() - 1;
    r2 = x*x + y*y;
  } while(r2 >= 1 || r2 == 0);
  double fac = sqrt(-2*log(r2)/r2);
  vector3d out(x*fac, y*fac, 0);
  do {
    x = 2*ran() - 1;
    y = 2*ran() - 1;
    r2 = x*x + y*y;
  } while(r2 >= 1 || r2 == 0);
  fac = sqrt(-2*log(r2)/r2);
  out[2]=x*fac;
  return out;
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

inline void print_all(const polyhedron *p) {
  if (debug) {
    for (int i=0; i<N; i++) {
      printf("%s %4i: (%6.2f, %6.2f, %6.2f)  [%5.2f, (%5.2f, %5.2f, %5.2f)] R: %4.2f\n", p[i].mypoly->name, i, p[i].pos[0], p[i].pos[1], p[i].pos[2], p[i].rot[0], p[i].rot[1], p[i].rot[2], p[i].rot[3], p[i].R);
      if (!in_cell(p[i]))
        printf("\t  Outside cell!\n");
      if (fabs(p[i].rot.normsquared() - 1.0) > 10e-15)
        printf("\t  Quaternion off! Normsquared - 1: %g\n", (p[i].rot.normsquared()-1.0));
      for (int j=i+1; j<N; j++)
        if (overlap(p[i], p[j]))
          printf("\t  Overlaps with %i!!!\n", j);
    }
  }
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
        printf("%s %4i: (%6.2f, %6.2f, %6.2f)  [%5.2f, (%5.2f, %5.2f, %5.2f)] R: %4.2f\n", p[i].mypoly->name, i, p[i].pos[0], p[i].pos[1], p[i].pos[2], p[i].rot[0], p[i].rot[1], p[i].rot[2], p[i].rot[3], p[i].R);
        if (!incell)
          printf("\t  Outside cell!\n");
        for (int j=0; j<N; j++) {
          if (j != i && overlap(p[i], p[j])) {
            printf("\t  Overlaps with %i", j);
            printf(": (%6.2f, %6.2f, %6.2f)  [%5.2f, (%5.2f, %5.2f, %5.2f)]\n", p[j].pos[0], p[j].pos[1], p[j].pos[2], p[j].rot[0], p[j].rot[1], p[j].rot[2], p[j].rot[3]);
          }
        }
      }
    }
  }
  fflush(stdout);
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

void define_shapes() {
  // Shapes, unrotated, at the origin, and circumscribed by a sphere of radius 1.0
  cube.nvertices = 8;
  cube.nfaces = 3;
  sprintf(cube.name, "cube");
  const double v_cube = 1.0/sqrt(3.0);
  cube.vertices[0] = vector3d(v_cube, v_cube, v_cube);
  cube.vertices[1] = vector3d(-v_cube, v_cube, v_cube);
  cube.vertices[2] = vector3d(-v_cube, -v_cube, v_cube);
  cube.vertices[3] = vector3d(v_cube, -v_cube, v_cube);
  cube.vertices[4] = vector3d(v_cube, v_cube, -v_cube);
  cube.vertices[5] = vector3d(-v_cube, v_cube, -v_cube);
  cube.vertices[6] = vector3d(-v_cube, -v_cube, -v_cube);
  cube.vertices[7] = vector3d(v_cube, -v_cube, -v_cube);

  cube.faces[0] = vector3d(1.0, 0, 0);
  cube.faces[1] = vector3d(0, 1.0, 0);
  cube.faces[2] = vector3d(0, 0, 1.0);



  tetrahedron.nvertices = 4;
  tetrahedron.nfaces = 4;
  sprintf(tetrahedron.name, "tetrahedron");
  tetrahedron.vertices[0] = vector3d(sqrt(2.0/3.0), 0, -sqrt(3.0));
  tetrahedron.vertices[1] = vector3d(-sqrt(2.0/3.0), 0, -sqrt(3.0));
  tetrahedron.vertices[2] = vector3d(0, sqrt(2.0/3.0), sqrt(3.0));
  tetrahedron.vertices[3] = vector3d(0, -sqrt(2.0/3.0), sqrt(3.0));

  tetrahedron.faces[0] = vector3d();
  tetrahedron.faces[1] = vector3d();
  tetrahedron.faces[2] = vector3d();
  tetrahedron.faces[3] = vector3d();




  truncated_tetrahedron.nvertices = 18;
  truncated_tetrahedron.nfaces = 4;




}
