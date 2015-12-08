#include "vector3d.h"

Rand random::my_rand = Rand(0);
unsigned long random::seedval = 0;

vector3d vector3d::ran(double scale) {
  double x, y, r2;
  do {
    x = 2*random::ran() - 1;
    y = 2*random::ran() - 1;
    r2 = x*x + y*y;
  } while(r2 >= 1 || r2 == 0);
  double fac = scale*sqrt(-2*log(r2)/r2);
  vector3d out(x*fac, y*fac, 0);
  static double extra_z = 0;
  if (extra_z) {
    // We have a z value left over from last time!
    out[2] = extra_z;
    extra_z = 0;
  } else {
    do {
      x = 2*random::ran() - 1;
      y = 2*random::ran() - 1;
      r2 = x*x + y*y;
    } while(r2 >= 1 || r2 == 0);
    fac = scale*sqrt(-2*log(r2)/r2);
    extra_z = y*fac; // Save this one for later!
    out[2]=x*fac;
  }
  return out;
}

vector3d vector3d::expran() {
  double x, y, z, r2;
  do {
    x = 2*random::ran() - 1;
    y = 2*random::ran() - 1;
    z = 2*random::ran() - 1;
    r2 = x*x + y*y + z*z;
  } while(r2 >= 1 || r2 == 0);
  double invdistance = -sqrt(r2)/log(random::ran());
  return vector3d(x*invdistance, y*invdistance, z*invdistance);
}

