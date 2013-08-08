#include <stdio.h>
//#include <time.h>
//#include "Monte-Carlo/monte-carlo.h"
#include "MersenneTwister.h"
#include "vector3d.h"
//#include <cassert>
#include <math.h>
//#include <stdlib.h>
//#include <string.h>
//#include <vector>
//using std::vector;

struct shape {
  vector3d pos;
  quaternion rot;
  double R;
};

// Global Constants

// Global "Constants" -- set at runtime then unchanged
int N = 10;
bool periodic[3] = {false, false, false};
bool walls[3] = {false, false, false};
double lenx, leny, lenz = 20;


// Shapes, unrotated, at the origin, and circumscribed by a sphere of radius 1.0
const vector3d tetrahedron[4] = {vector3d(sqrt(2.0/3.0), 0, -sqrt(3.0)),
                                 vector3d(-sqrt(2.0/3.0), 0, -sqrt(3.0)),
                                 vector3d(0, sqrt(2.0/3.0), sqrt(3.0)),
                                 vector3d(0, -sqrt(2.0/3.0), sqrt(3.0))};
const double v_cube = 1.0/sqrt(3.0);
const vector3d cube[8] = {vector3d(v_cube, v_cube, v_cube),
                          vector3d(-v_cube, v_cube, v_cube),
                          vector3d(-v_cube, -v_cube, v_cube),
                          vector3d(v_cube, -v_cube, v_cube),
                          vector3d(v_cube, v_cube, -v_cube),
                          vector3d(-v_cube, v_cube, -v_cube),
                          vector3d(-v_cube, -v_cube, -v_cube),
                          vector3d(v_cube, -v_cube, -v_cube)};


// Check whether two polyhedra overlap
bool overlap(const shape &a, const shape &b);

// Return the vector pointing from a to b, accounting for periodic
// boundaries
vector3d periodic_diff(const vector3d &a, const vector3d  &b);

// Move v in a random direction by a distance determined by a gaussian distribution
vector3d move(const vector3d &v, double scale);

// Rotate q in a random direction by an amount determined by a gaussian distribution
quaternion rotate(const quaternion &q, double scale);

// Rotate the vector v about the axis of vector part of q
// by an amount equal to the scalar part of q
vector3d rotate_vector(const vector3d &v, const quaternion &q);

// Generate a random number in the range [0, 1) using a fixed seed
double ran();

// Generate a random point in a gaussian distribution with sigma = 1.0
vector3d ran3();

// Generate a random quaternion, with vector part normalized
// and scalar part in a gaussian distribution with sigma = 1.0
quaternion ran4();

// If v is outside the cell, and there are periodic boundary condition(s), it is
// moved (appropriately) into the cell
inline vector3d fix_periodic(vector3d newv);



int main() {
  shape *shapes = new shape[N];
  for (int i=0; i<N; i++) {
    shapes[i].rot[3] = 1.0;
    shapes[i].rot[0] = M_PI/2.0*i;
    shapes[i].pos[0] = 1.0;
    shapes[i].pos = rotate_vector(shapes[i].pos, shapes[i].rot);
    printf("(%.2f, %.2f, %.2g)\t[%.2g, (%.2g, %.2g, %.2g)]\n", shapes[i].pos[0], shapes[i].pos[1], shapes[i].pos[2], shapes[i].rot[0], shapes[i].rot[1], shapes[i].rot[2], shapes[i].rot[3]);
  }
  return 0;
}


bool overlap(const shape &a, const shape &b) {
  const vector3d ab = periodic_diff(a.pos, b.pos);
  if (ab.norm() > a.R + b.R)
    return false;
  return true;
}

vector3d periodic_diff(const vector3d &a, const vector3d &b) {
  vector3d v = b - a;
  if (periodic[0]) {
    while (v.x() > lenx/2.0)
      v[0] -= lenx;
    while (v.x() < -lenx/2.0)
      v[0] += lenx;
  }
  if (periodic[1]) {
    while (v.y() > leny/2.0)
      v[1] -= leny;
    while (v.y() < -leny/2.0)
      v[1] += leny;
  }
  if (periodic[2]) {
    while (v.z() > lenz/2.0)
      v[2] -= lenz;
    while (v.z() < -lenz/2.0)
      v[2] += lenz;
  }
  return v;
}

vector3d move(const vector3d &v, double scale) {
  const vector3d newv = v+scale*ran3();
  return fix_periodic(newv);
}

quaternion rotate(const quaternion &q, double scale) {
  quaternion q2 = ran4();
  q2[0] *= scale;
  return q2*q;
}

vector3d rotate_vector(const vector3d &v, const quaternion &q) {
  const quaternion p(0, v[0], v[1], v[2]);
  const quaternion product = q*p*q.conj();
  return vector3d(product[1], product[2], product[3]);
}


double ran() {
  const long unsigned int x = 0;
  static MTRand my_mtrand(x); // always use the same random number generator (for debugging)!
  return my_mtrand.randExc(); // which is the range of [0,1)
}

vector3d ran3() {
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

quaternion ran4() {
  double x, y, r2;
  do {
    x = 2*ran() - 1;
    y = 2*ran() - 1;
    r2 = x*x + y*y;
  } while(r2 >= 1 || r2 == 0);
  double fac = sqrt(-2*log(r2)/r2);
  quaternion out(x*fac, y*fac, 0, 0);
  do {
    x = 2*ran() - 1;
    y = 2*ran() - 1;
    r2 = x*x + y*y;
  } while(r2 >= 1 || r2 == 0);
  fac = sqrt(-2*log(r2)/r2);
  out[2] = x*fac;
  out[3] = y*fac;
  return out.normalize_vector();
}

inline vector3d fix_periodic(vector3d newv) {
  if (periodic[0] || walls[0]) {
    while (newv[0] > lenx/2)
      newv[0] -= lenx;
    while (newv[0] < -lenx/2)
      newv[0] += lenx;
  }
  if (periodic[1] || walls[1]) {
    while (newv[1] > leny/2)
      newv[1] -= leny;
    while (newv[1] < -leny/2)
      newv[1] += leny;
  }
  if (periodic[2] || walls[2]) {
    while (newv[2] > lenz/2)
      newv[2] -= lenz;
    while (newv[2] < -lenz/2)
      newv[2] += lenz;
  }
  return newv;
}
