#include <math.h>

#pragma once

class vector3d {
 public:
  vector3d() {
    coords[0] = coords[1] = coords[2] = 0; }
  vector3d(const double x, const double y, const double z) {
    coords[0] = x; coords[1] = y; coords[2] = z; }
  vector3d(const vector3d &v) {
    coords[0] = v.coords[0]; coords[1] = v.coords[1]; coords[2] = v.coords[2];}

  double x() const { return coords[0]; }
  double y() const { return coords[1]; }
  double z() const { return coords[2]; }

  double x(const double newx) { coords[0] = newx; return coords[0]; }
  double y(const double newy) { coords[1] = newy; return coords[1]; }
  double z(const double newz) { coords[2] = newz; return coords[2]; }

  vector3d operator=(const vector3d &v) {
    coords[0] = v.coords[0]; coords[1] = v.coords[1]; coords[2] = v.coords[2];
    return *this; }

  vector3d operator+(const vector3d &v) const {
    return vector3d(coords[0]+v.coords[0], coords[1]+v.coords[1], coords[2]+v.coords[2]); }
  vector3d operator+=(const vector3d &v) {
    coords[0] += v.coords[0]; coords[1] += v.coords[1]; coords[2] += v.coords[2];
    return *this; }

  vector3d operator-(const vector3d &v) const {
    return vector3d(coords[0]-v.coords[0], coords[1]-v.coords[1], coords[2]-v.coords[2]); }
  vector3d operator-=(const vector3d &v) {
    coords[0] -= v.coords[0]; coords[1] -= v.coords[1]; coords[2] -= v.coords[2];
    return *this; }

  vector3d operator*(const double scalar) const {
    return vector3d(scalar*coords[0], scalar*coords[1], scalar*coords[2]); }
  vector3d operator*=(const double scalar) {
    coords[0] *= scalar; coords[1] *= scalar; coords[2] *= scalar; return *this; }

  vector3d operator/(const double scalar) const {
    return vector3d(coords[0]/scalar, coords[1]/scalar, coords[2]/scalar); }
  vector3d operator/=(const double scalar) {
    coords[0] /= scalar; coords[1] /= scalar; coords[2] /= scalar; return *this; }

  bool operator ==(const vector3d &v) const {
    return ((coords[0] == v.coords[0]) && (coords[1] == v.coords[1]) &&
            (coords[2] == v.coords[2])); }
  bool operator !=(const vector3d &v) const {
    return !(*this == v); }

  double &operator[](const unsigned int i) { return coords[i]; }
  const double operator[](const unsigned int i) const { return coords[i]; }

  double dot(const vector3d &v) const {
    return coords[0]*v.coords[0] + coords[1]*v.coords[1] + coords[2]*v.coords[2]; }
  vector3d cross(const vector3d &v) const {
    return vector3d(coords[1]*v.coords[2] - coords[2]*v.coords[1],
                   coords[0]*v.coords[2] - coords[2]*v.coords[0],
                   coords[0]*v.coords[1] - coords[1]*v.coords[2]); }

  double norm() const {
    return sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]); }
  double normsquared() const {
    return coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]; }
  vector3d normalized() const {
    return *this/this->norm(); }

  void tostr(char str[]) const {
    sprintf(str, "(%6.2f, %6.2f, %6.2f)", coords[0], coords[1], coords[2]);
  }

 private:
  double coords[3];
};

const vector3d operator*(const double scalar, const vector3d &v) {
  return v*scalar; }


class quaternion {
 public:
  quaternion() {
    coords[0] = 1.0; coords[1] = coords[2] = coords[3] = 0; }
  quaternion(const double w, const double x, const double y, const double z) {
    coords[0] = w; coords[1] = x; coords[2] = y; coords[3] = z; }
  quaternion(const quaternion &q) {
    coords[0] = q.coords[0]; coords[1] = q.coords[1];
    coords[2] = q.coords[2]; coords[3] = q.coords[3]; }
  quaternion(const double w, const vector3d &v) {
    coords[0] = w; coords[1] = v[0]; coords[2] = v[1]; coords[3] = v[2]; }

  double w() const { return coords[0]; }
  double x() const { return coords[1]; }
  double y() const { return coords[2]; }
  double z() const { return coords[3]; }

  double w(const double neww) { coords[0] = neww; return coords[0]; }
  double x(const double newx) { coords[1] = newx; return coords[1]; }
  double y(const double newy) { coords[2] = newy; return coords[2]; }
  double z(const double newz) { coords[3] = newz; return coords[3]; }

  quaternion operator=(const quaternion &q) {
    coords[0] = q.coords[0]; coords[1] = q.coords[1];
    coords[2] = q.coords[2]; coords[3] = q.coords[3];
    return *this; }

  quaternion operator+(const quaternion &q) const {
    return quaternion(coords[0]+q.coords[0], coords[1]+q.coords[1],
                      coords[2]+q.coords[2], coords[3]+q.coords[3]); }
  quaternion operator+=(const quaternion &q) {
    coords[0] += q.coords[0]; coords[1] += q.coords[1];
    coords[2] += q.coords[2]; coords[3] += q.coords[3];
    return *this; }

  quaternion operator-(const quaternion &q) const {
    return quaternion(coords[0]-q.coords[0], coords[1]-q.coords[1],
                      coords[2]-q.coords[2], coords[3]-q.coords[3]); }
  quaternion operator-=(const quaternion &q) {
    coords[0] -= q.coords[0]; coords[1] -= q.coords[1];
    coords[2] -= q.coords[2]; coords[3] -= q.coords[3];
    return *this; }

  quaternion operator*(const double scalar) const {
    return quaternion(scalar*coords[0], scalar*coords[1],
                      scalar*coords[2], scalar*coords[3]); }
  quaternion operator*=(const double scalar) {
    coords[0] *= scalar; coords[1] *= scalar; coords[2] *= scalar; coords[3] *= scalar;
    return *this; }

  quaternion operator*(const quaternion &q) const {
    return quaternion(coords[0]*q.coords[0] - coords[1]*q.coords[1] -
                      coords[2]*q.coords[2] - coords[3]*q.coords[3],
                      coords[0]*q.coords[1] + coords[1]*q.coords[0] +
                      coords[2]*q.coords[3] - coords[3]*q.coords[2],
                      coords[0]*q.coords[2] - coords[1]*q.coords[3] +
                      coords[2]*q.coords[0] + coords[3]*q.coords[1],
                      coords[0]*q.coords[3] + coords[1]*q.coords[2] -
                      coords[2]*q.coords[1] + coords[3]*q.coords[0]); }
  quaternion operator*=(const quaternion &q) {
    *this = (*this)*q; return *this; }

  quaternion operator/(const double scalar) const {
    return quaternion(coords[0]/scalar, coords[1]/scalar,
                      coords[2]/scalar, coords[3]/scalar); }
  quaternion operator/=(const double scalar) {
    coords[0]/=scalar; coords[1]/=scalar; coords[2]/=scalar; coords[3]/=scalar;
    return *this; }

  bool operator ==(const quaternion &q) const {
    return ((coords[0] == q.coords[0]) && (coords[1] == q.coords[1]) &&
            (coords[2] == q.coords[2]) && (coords[3] == q.coords[3])); }
  bool operator !=(const quaternion &q) const {
    return !(*this == q); }

  double &operator[](const unsigned int i) { return coords[i]; }
  const double &operator[](const unsigned int i) const { return coords[i]; }

  quaternion conj() const {
    return quaternion(coords[0], -coords[1], -coords[2], -coords[3]); }

  double norm() const {
    return sqrt(coords[0]*coords[0] + coords[1]*coords[1] +
                coords[2]*coords[2] + coords[3]*coords[3]); }
  double normsquared() const {
    return coords[0]*coords[0] + coords[1]*coords[1] +
                coords[2]*coords[2] + coords[3]*coords[3]; }
  quaternion normalized_vector() const {
    const double norm = sqrt(coords[1]*coords[1] + coords[2]*coords[2] +
                             coords[3]*coords[3]);
    return quaternion(coords[0], coords[1]/norm, coords[2]/norm, coords[3]/norm); }
  quaternion normalized() const {
    return *this/this->norm();
  }

  void tostr(char str[]) const {
    sprintf(str, "[%6.2f, (%6.2f, %6.2f, %6.2f)]", coords[0], coords[1], coords[2], coords[3]); }

 private:
  double coords[4]; // coords[0] = scalar part, rest is vector part
};

const quaternion operator*(const double scalar, const quaternion &q) {
  return q*scalar; }
