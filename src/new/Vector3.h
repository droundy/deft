// -*- mode: C++; -*-

#pragma once

#include <cassert>

class Vector3 {
public:
  Vector3(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }
  Vector3(const Vector3 &a) {
    for (int i=0; i<3; i++) {
      data[i] = a.data[i];
    }
  }
  void operator=(const Vector3 &a) {
    for (int i=0; i<3; i++) {
      data[i] = a.data[i];
    }
  }

  // operator[] is a pair of "safe" operator for indexing a vector
  double operator[](int i) const {
    assert(i >= 0 && i < 3);
    return data[i];
  };
  double &operator[](int i) {
    assert(i >= 0 && i < 3);
    return data[i];
  };

  // Below are the arithmetic operators.  These should be made fast.
  Vector3 operator+(const Vector3 &a) const {
    Vector3 out;
    for (int i=0; i<3; i++) {
      out[i] = data[i] + a.data[i];
    }
    return out;
  }
  Vector3 operator-(const Vector3 &a) const {
    Vector3 out;
    for (int i=0; i<3; i++) {
      out[i] = data[i] - a.data[i];
    }
    return out;
  }
  Vector3 operator-(const Vector3 &a) const {
    Vector3 out;
    for (int i=0; i<3; i++) {
      out[i] = data[i] - a.data[i];
    }
    return out;
  }
  Vector3 operator*(double a) const {
    Vector3 out;
    for (int i=0; i<3; i++) {
      out[i] = a*data[i];
    }
    return out;
  }
  friend Vector3 operator*(double a, const Vector3 &b);

  void operator+=(const Vector3 &a) {
    for (int i=0; i<3; i++) {
      data[i] += a.data[i];
    }
    return out;
  }
  void operator-=(const Vector3 &a) {
    for (int i=0; i<3; i++) {
      data[i] -= a.data[i];
    }
    return out;
  }
  void operator*=(double a) {
    for (int i=0; i<3; i++) {
      data[i] *= a;
    }
  }

  double dot(const Vector3 &b) const {
    return data[0]*b.data[0] + data[1]*b.data[1] + data[2]*b.data[2];
  }
  Vector3 cross(const Vector3 &b) const {
    Vector3 out;
    for (int i=0; i<3; i++) {
      out.data[i] = data[(i+1)%3]*b.data[(i+2)%3] - data[(i+2)%3]*b.data[(i+1)%3];
    }
    return out;
  }
private:
  double data[3];
};

inline Vector3 operator*(double a, const Vector3 &b) {
  Vector3 out;
  for (int i=0; i<3; i++) {
    out.data[i] = a*b.data[i];
  }
  return out;
}
