#pragma once

#include <eigen2/Eigen/Core>
#include <eigen2/Eigen/LU>
#include <assert.h>
#include <math.h>

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class Relative : public Vector3d {
public:
  explicit Relative(const Vector3d &x) : Vector3d(x) {}
  Relative(double x,double y,double z) : Vector3d(x,y,z) {}
};

class Cartesian : public Vector3d {
public:
  explicit Cartesian(const Vector3d &x) : Vector3d(x) {}
  Cartesian(double x,double y,double z) : Vector3d(x,y,z) {}
};

class Reciprocal : public Vector3d {
public:
  explicit Reciprocal(const Vector3d &x) : Vector3d(x) {}
  Reciprocal(double x,double y,double z) : Vector3d(x,y,z) {}
};

class RelativeReciprocal : public Vector3d {
public:
  explicit RelativeReciprocal(const Vector3d &x) : Vector3d(x) {}
  RelativeReciprocal(double x,double y,double z) : Vector3d(x,y,z) {}
};

class Lattice {
  Matrix3d R, G;
  double vol;
public:
  Lattice(Cartesian a1, Cartesian a2, Cartesian a3) {
    R.col(0) = a1; R.col(1) = a2; R.col(2) = a3;
    G = R.inverse();
    vol = R.determinant();
    assert(volume() > 0);
  };
  Cartesian toCartesian(Relative x) const {
    return Cartesian(R*x);
  };
  Relative toRelative(Cartesian x) const {
    return Relative(G*x);
  };
  Reciprocal toReciprocal(RelativeReciprocal x) const {
    return Reciprocal(G.transpose()*x);
  };
  RelativeReciprocal toRelativeReciprocal(Reciprocal x) const {
    return RelativeReciprocal(R.transpose()*x);
  };
  double volume() const { return vol; };
};

// Define the dot products between different types that make sense.
double operator*(Cartesian r, Reciprocal k) {
  return r.dot(k);
}
double operator*(Cartesian r, Cartesian r2) {
  return r.dot(r2);
}
double operator*(Reciprocal k, Cartesian r) {
  return r.dot(k);
}
double operator*(Reciprocal k, Reciprocal k2) {
  return k.dot(k2);
}
double operator*(Relative r, RelativeReciprocal k) {
  return (2*M_PI)*r.dot(k);
}
double operator*(RelativeReciprocal k, Relative r) {
  return r.dot(k);
}
