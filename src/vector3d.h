#include <math.h>
#include "MersenneTwister.h"

#pragma once

struct random {
  static void seed(unsigned long seedval) { my_mtrand = MTRand(seedval); }
  static double ran() {return my_mtrand.randExc(); }
private:
  static MTRand my_mtrand;
};

class vector3d {
 public:
  double x;
  double y;
  double z;

  vector3d() {
    x = y = z = 0; }
  vector3d(const double newx, const double newy, const double newz) {
    x = newx; y = newy; z = newz; }
  vector3d(const vector3d &v) {
    x = v.x; y = v.y; z = v.z;}

  vector3d operator=(const vector3d &v) {
    x = v.x; y = v.y; z = v.z;
    return *this; }

  vector3d operator+(const vector3d &v) const {
    return vector3d(x+v.x, y+v.y, z+v.z); }
  vector3d operator+=(const vector3d &v) {
    x += v.x; y += v.y; z += v.z;
    return *this; }

  vector3d operator-(const vector3d &v) const {
    return vector3d(x-v.x, y-v.y, z-v.z); }
  vector3d operator-=(const vector3d &v) {
    x -= v.x; y -= v.y; z -= v.z;
    return *this; }

  vector3d operator*(const double scalar) const {
    return vector3d(scalar*x, scalar*y, scalar*z); }
  vector3d operator*=(const double scalar) {
    x *= scalar; y *= scalar; z *= scalar;
    return *this; }

  vector3d operator/(const double scalar) const {
    return vector3d(x/scalar, y/scalar, z/scalar); }
  vector3d operator/=(const double scalar) {
    x /= scalar; y /= scalar; z /= scalar;
    return *this; }

  bool operator ==(const vector3d &v) const {
    return ((x == v.x) && (y == v.y) &&
            (z == v.z)); }
  bool operator !=(const vector3d &v) const {
    return !(*this == v); }

  double &operator[](const unsigned int i) {
    switch(i) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    }
  }
  const double operator[](const unsigned int i) const {
    switch(i) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    }
  }

  double dot(const vector3d &v) const {
    return x*v.x + y*v.y + z*v.z; }
  vector3d cross(const vector3d &v) const {
    return vector3d(y*v.z - z*v.y, x*v.z - z*v.x, x*v.y - y*v.z); }

  double norm() const {
    return sqrt(x*x + y*y + z*z); }
  double normsquared() const {
    return x*x + y*y + z*z; }
  vector3d normalized() const {
    return *this/this->norm(); }

  void tostr(char str[]) const {
    sprintf(str, "(%6.2f, %6.2f, %6.2f)", x, y, z);
  }

  static vector3d ran(double scale) {
    double x, y, r2;
    do {
      x = 2*random::ran() - 1;
      y = 2*random::ran() - 1;
      r2 = x*x + y*y;
    } while(r2 >= 1 || r2 == 0);
    double fac = sqrt(-2*log(r2)/r2);
    vector3d out(x*fac, y*fac, 0);
    do {
      x = 2*random::ran() - 1;
      y = 2*random::ran() - 1;
      r2 = x*x + y*y;
    } while(r2 >= 1 || r2 == 0);
    fac = sqrt(-2*log(r2)/r2);
    out[2]=x*fac;
    return out;
  }
};

const vector3d operator*(const double scalar, const vector3d &v) {
  return v*scalar; }

class rotation {
 public:
  double w;
  double x;
  double y;
  double z;

  rotation() { w = 1.0; x = y = z = 0; }
  rotation(const rotation &q) { w = q.w; x = q.x; y = q.y; z = q.z; }

  rotation operator=(const rotation &q) {
    w = q.w; x = q.x; y = q.y; z = q.z;
    return *this; }

  rotation operator*(const rotation &q) const {
    return rotation(w*q.w - x*q.x - y*q.y - z*q.z,
                    w*q.x + x*q.w + y*q.z - z*q.y,
                    w*q.y - x*q.z + y*q.w + z*q.x,
                    w*q.z + x*q.y - y*q.x + z*q.w).normalized(); }
  rotation operator*=(const rotation &q) {
    *this = (*this)*q; return *this; }

  bool operator ==(const rotation &q) const {
    return ((w == q.w) && (x == q.x) &&
            (y == q.y) && (z == q.z)); }
  bool operator !=(const rotation &q) const {
    return !(*this == q); }

  vector3d rotate_vector(const vector3d &v) {
    const rotation product = (*this)*rotation(0, v)*this->conj();
    return vector3d(product.x, product.y, product.z); }

  void tostr(char str[]) const {
    const double theta = 2*acos(w);
    const double fac = 1.0/sin(theta/2.0);
    sprintf(str, "[%6.2f, (%6.2f, %6.2f, %6.2f)]", theta, x*fac, y*fac, z*fac); }

  static rotation ran() {
    double x, y, z, r2;
    const double theta = 2.0*M_PI*random::ran();
    do {
      x = 2*random::ran() - 1;
      y = 2*random::ran() - 1;
      z = 2*random::ran() - 1;
      r2 = x*x + y*y + z*z;
    } while(r2 >= 1 || r2 == 0);
    const rotation rot(cos(theta/2), vector3d(x, y, z)*sin(theta/2)/sqrt(r2));
    return rot.normalized();
  }

  static rotation ran(double angwidth) {
    double x, y, z, r2, sintheta_over2;
    do {
      do {
        x = 2*random::ran() - 1;
        y = 2*random::ran() - 1;
        r2 = x*x + y*y;
      } while (r2 >= 1 || r2 == 0);
    const double fac = sqrt(-2*log(r2)/r2);
    sintheta_over2 = fac*x*angwidth;
    } while (sintheta_over2 <= -1 or sintheta_over2 >= 1);
    const double costheta_over2 = sqrt(1 - sintheta_over2*sintheta_over2);
    do {
      x = 2*random::ran() - 1;
      y = 2*random::ran() - 1;
      z = 2*random::ran() - 1;
      r2 = x*x + y*y + z*z;
    } while(r2 >= 1 || r2 == 0);
    const double vfac = sintheta_over2/sqrt(r2);
    const rotation rot(costheta_over2, x*vfac, y*vfac, z*vfac);
    return rot.normalized();
  }

 private:
  rotation(const double neww, const double newx, const double newy, const double newz) {
    w = neww; x = newx; y = newy; z = newz; }
  rotation(const double neww, const vector3d &v) {
    w = neww; x = v.x; y = v.y; z = v.z; }

  rotation operator/(const double scalar) const {
    return rotation(w/scalar, x/scalar, y/scalar, z/scalar); }
  rotation operator/=(const double scalar) {
    w/=scalar; x/=scalar; y/=scalar; z/=scalar;
    return *this; }

  rotation conj() const {
    return rotation(w, -x, -y, -z); }

  rotation normalized() const {
    return *this/sqrt(w*w + x*x + y*y + z*z);
  }
};

class quaternion {
 public:
  double w;
  double x;
  double y;
  double z;

  quaternion() { w = 1.0; x = y = z = 0; }
  quaternion(const double neww, const double newx, const double newy, const double newz) {
    w = neww; x = newx; y = newy; z = newz; }
  quaternion(const quaternion &q) {
    w = q.w; x = q.x; y = q.y; z = q.z; }
  quaternion(const double neww, const vector3d &v) {
    w = neww; x = v.x; y = v.y; z = v.z; }

  quaternion operator=(const quaternion &q) {
    w = q.w; x = q.x; y = q.y; z = q.z;
    return *this; }

  quaternion operator+(const quaternion &q) const {
    return quaternion(w+q.w, x+q.x, y+q.y, z+q.z); }
  quaternion operator+=(const quaternion &q) {
    w += q.w; x += q.x; y += q.y; z += q.z;
    return *this; }

  quaternion operator-(const quaternion &q) const {
    return quaternion(w-q.w, x-q.x, y-q.y, z-q.z); }
  quaternion operator-=(const quaternion &q) {
    w -= q.w; x -= q.x; y -= q.y; z -= q.z;
    return *this; }

  quaternion operator*(const double scalar) const {
    return quaternion(scalar*w, scalar*x,
                      scalar*y, scalar*z); }
  quaternion operator*=(const double scalar) {
    w *= scalar; x *= scalar; y *= scalar; z *= scalar;
    return *this; }

  quaternion operator*(const quaternion &q) const {
    return quaternion(w*q.w - x*q.x - y*q.y - z*q.z,
                      w*q.x + x*q.w + y*q.z - z*q.y,
                      w*q.y - x*q.z + y*q.w + z*q.x,
                      w*q.z + x*q.y - y*q.x + z*q.w); }
  quaternion operator*=(const quaternion &q) {
    *this = (*this)*q; return *this; }

  quaternion operator/(const double scalar) const {
    return quaternion(w/scalar, x/scalar,
                      y/scalar, z/scalar); }
  quaternion operator/=(const double scalar) {
    w/=scalar; x/=scalar; y/=scalar; z/=scalar;
    return *this; }

  bool operator ==(const quaternion &q) const {
    return ((w == q.w) && (x == q.x) &&
            (y == q.y) && (z == q.z)); }
  bool operator !=(const quaternion &q) const {
    return !(*this == q); }

  double &operator[](const unsigned int i) {
    switch(i) {
    case 0: return w;
    case 1: return x;
    case 2: return y;
    case 3: return z;
    }
  }
  const double operator[](const unsigned int i) const {
    switch(i) {
    case 0: return w;
    case 1: return x;
    case 2: return y;
    case 3: return z;
    }
  }

  quaternion conj() const {
    return quaternion(w, -x, -y, -z); }

  double norm() const {
    return sqrt(w*w + x*x + y*y + z*z); }
  double normsquared() const {
    return w*w + x*x + y*y + z*z; }
  quaternion normalized_vector() const {
    const double norm = sqrt(x*x + y*y + z*z);
    return quaternion(w, x/norm, y/norm, z/norm); }
  quaternion normalized() const {
    return *this/this->norm();
  }

  void tostr(char str[]) const {
    sprintf(str, "[%6.2f, (%6.2f, %6.2f, %6.2f)]", w, x, y, z); }
};

const quaternion operator*(const double scalar, const quaternion &q) {
  return q*scalar; }
