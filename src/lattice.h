#pragma once

#include <eigen2/Eigen/Core>
#include <eigen2/Eigen/LU>
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
  Matrix3d R, G;
  double vol;
public:
  Lattice(Cartesian a1, Cartesian a2, Cartesian a3) {
    R.col(0) = a1; R.col(1) = a2; R.col(2) = a3;
    G = R.inverse();
    vol = R.determinant();
    assert(volume() > 0);
  };
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
  void reorientBasis(Cartesian zdir) {
    zdir = zdir / zdir.norm();
    Cartesian a1new(1e10, 0, 0.3e9),
              a2new(1.3e9, 1e10, 0.7e9),
              a3new(1.23e9, 1e9, 0.73e9);
    bool got1 = false, got2 = false;
    for (int i=-3; i<=3; i++) {
      for (int j=-3; j<=3; j++) {
        for (int k=-3; k<=3; k++) {
          Cartesian here(i*R.col(0) + j*R.col(1) + k*R.col(2));
          double sin1 = here.cross(a1new).norm()/(here.norm()*a1new.norm());
          double sin2 = here.cross(a2new).norm()/(here.norm()*a2new.norm());
          double cosz = here.dot(zdir)/here.norm();
          if (fabs(cosz) <= 1e-7) {
            if (!got1) {
              a1new = here;
              got1 = true;
            } else if (!got2 && sin1 >= 1e-7) {
              a2new = here;
              got2 = true;
            } else if (here.norm() < a1new.norm() && sin2 >= 1e-7) {
              a1new = here;
            } else if (here.norm() < a2new.norm() && sin1 >= 1e-7) {
              a2new = here;
            }
          } else {
            if (here.norm() > 0 &&
                fabs(here.dot(zdir)) < fabs(a3new.dot(zdir)))
              a3new = here;
          }
        }
      }
    }
    // set up the new lattice...
    double vol0 = volume();
    R.col(0) = a1new; R.col(1) = a2new; R.col(2) = a3new;
    if (R.determinant() < 0) R.col(2) = -a3new;
    G = R.inverse();
    vol = R.determinant();
    assert(volume() > 0);
    assert(fabs(volume()-vol0)/vol0 < 1e-7);
  }
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
