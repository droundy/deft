// -*- mode: C++; -*-

#pragma once

#pragma GCC diagnostic push
#if __GCC__ > 5
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Core>
#pragma GCC diagnostic pop

#include <assert.h>
#include <math.h>

// import most common Eigen types 
USING_PART_OF_NAMESPACE_EIGEN

class Relative : public Vector3d {
public:
  Relative(void) : Vector3d() {}
  explicit Relative(const Vector3d &x) : Vector3d(x) {}
  Relative(double x,double y,double z) : Vector3d(x,y,z) {}

  // We need to define this for our object to work
  typedef Eigen::Vector3d Base;
  template<typename OtherDerived>
  Relative &operator= (const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }
};

class Cartesian : public Vector3d {
public:
  Cartesian(void) : Vector3d() {}
  explicit Cartesian(const Vector3d &x) : Vector3d(x) {}
  Cartesian(double x,double y,double z) : Vector3d(x,y,z) {}

  // We need to define this for our object to work
  typedef Eigen::Vector3d Base;
  template<typename OtherDerived>
  Cartesian &operator= (const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }
};

class Reciprocal : public Vector3d {
public:
  Reciprocal(void) : Vector3d() {}
  explicit Reciprocal(const Vector3d &x) : Vector3d(x) {}
  Reciprocal(double x,double y,double z) : Vector3d(x,y,z) {}

  // We need to define this for our object to work
  typedef Eigen::Vector3d Base;
  template<typename OtherDerived>
  Reciprocal &operator= (const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }
};

class RelativeReciprocal : public Vector3d {
public:
  RelativeReciprocal(void) : Vector3d() {}
  explicit RelativeReciprocal(const Vector3d &x) : Vector3d(x) {}
  RelativeReciprocal(double x,double y,double z) : Vector3d(x,y,z) {}

  // We need to define this for our object to work
  typedef Eigen::Vector3d Base;
  template<typename OtherDerived>
  RelativeReciprocal &operator= (const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }
};

class Lattice {
  Matrix3d R, Rinv, minimalR, minimalRinv;
  double vol;
public:
  Lattice(Cartesian a1, Cartesian a2, Cartesian a3);
  Cartesian a1() const { return Cartesian(R.col(0)); }
  Cartesian a2() const { return Cartesian(R.col(1)); }
  Cartesian a3() const { return Cartesian(R.col(2)); }
  Cartesian toCartesian(Relative x) const {
    return Cartesian(R*x);
  };
  Cartesian round(Cartesian x) const {
    return toCartesian(round(toRelative(x)));
  };
  Relative round(Relative x) const {
    return Relative(floor(x(0)+0.5), floor(x(1)+0.5), floor(x(2)+0.5));
  };

  // reorientBasis modifies the basis so the first two vectors are
  // orthogonal to zdir.
  void reorientBasis(Cartesian zdir);
  Cartesian wignerSeitz(Cartesian) const;
  Reciprocal brillouinZone(Reciprocal) const;
  Relative toRelative(Cartesian x) const {
    return Relative(Rinv*x);
  };
  Reciprocal toReciprocal(RelativeReciprocal x) const {
    return Reciprocal(2*M_PI*Rinv.transpose()*x);
  };
  RelativeReciprocal toRelativeReciprocal(Reciprocal x) const {
    return RelativeReciprocal((1.0/2/M_PI)*R.transpose()*x);
  };
  double volume() const { return vol; };
};

// Define the dot products between different types that make sense.
inline double operator*(Cartesian r, Reciprocal k) {
  return r.dot(k);
}
inline double operator*(Cartesian r, Cartesian r2) {
  return r.dot(r2);
}
inline double operator*(Reciprocal k, Cartesian r) {
  return r.dot(k);
}
inline double operator*(Reciprocal k, Reciprocal k2) {
  return k.dot(k2);
}
inline double operator*(Relative r, RelativeReciprocal k) {
  return (2*M_PI)*r.dot(k);
}
inline double operator*(RelativeReciprocal k, Relative r) {
  return r.dot(k);
}
