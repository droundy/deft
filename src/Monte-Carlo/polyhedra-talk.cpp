#include <stdio.h>
#include <time.h>
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <popt.h>
#include "handymath.h"
#include "vector3d.h"
#include "Monte-Carlo/polyhedra.h"

// Saves the vertices of all polyhedra to a file.
inline void save_locations(const polyhedron *p, int N, const char *fname,
                           const double len[3], const char *comment="");


int main(int argc, const char *argv[]) {
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
