#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"
#include "vector3d.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
using std::vector;

struct shape {
  vector3d pos;
  quaternion rot;
};

// Global constants
const double R = 0; // the radius of the sphere that circumscribes each shape

// Global Variables - set at runtime then unchanged
int N = 10;
bool periodic[3] = {false, false, false};
double lenx, leny, lenz = 20;

// Shapes, unrotated and at the origin
const vector3d tetrahedron[4] = {vector3d(sqrt(2.0/3.0)*R, 0, -sqrt(3.0)),
                                vector3d(-sqrt(2.0/3.0)*R, 0, -sqrt(3.0)),
                                vector3d(0, sqrt(2.0/3.0)*R, sqrt(3.0)),
                                vector3d(0, -sqrt(2.0/3.0)*R, sqrt(3.0))};
// Functions
double overlap(const shape &, const shape &);
vector3d periodic_diff(const vector3d &, const vector3d  &);

int main() {
  shape *shapes = new shape[N];
  for (int i=0; i<N; i++) {
    shapes[i].pos[0] += i*1;
  }
  printf("%g\n", overlap(shapes[0], shapes[1]));
  return 0;
}

double overlap(const shape &a, const shape &b) {
  const vector3d ab = periodic_diff(a.pos, b.pos);
  if (ab.norm() > 1)
    return 0;
  return 1;
}

vector3d periodic_diff(const vector3d &a, const vector3d &b){
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
