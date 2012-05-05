// -*- mode: C++; -*-

#pragma once

#include <cassert>
#include <string.h>

class Vector {
public:
  Vector() : size(0), data(0) {}
  Vector(int sz) : size(sz), data(new double[size]) {}
  Vector(const Vector &a) : size(a.size), data(new double[size]) {
    memcpy(data, a.data, size*sizeof(double)); // faster than manual loop?
  }
  ~Vector() { delete[] data; }
  void operator=(const Vector &a) {
    if (size == 0) {
      size = a.size;
      data = new double[size];
    }
    assert(size = a.size);
    memcpy(data, a.data, size*sizeof(double)); // faster than manual loop?
  }

  // operator[] is a pair of "safe" operator for indexing a vector
  double operator[](int i) const {
    assert(i >= 0);
    assert(i < size);
    return data[i];
  };
  double &operator[](int i) {
    assert(i >= 0);
    assert(i < size);
    return data[i];
  };

  // Below are the arithmetic operators.  These should be made fast.
  Vector operator+(const Vector &a) const {
    assert(size == a.size);
    Vector out(size);
    for (int i=0; i<size; i++) {
      out[i] = data[i] + a.data[i];
    }
    return out;
  }
  Vector operator-(const Vector &a) const {
    assert(size == a.size);
    Vector out(size);
    for (int i=0; i<size; i++) {
      out[i] = data[i] - a.data[i];
    }
    return out;
  }
  Vector operator*(double a) const {
    Vector out(size);
    for (int i=0; i<size; i++) {
      out[i] = a*data[i];
    }
    return out;
  }
  friend Vector operator*(double a, const Vector &b);

  void operator+=(const Vector &a) {
    assert(size == a.size);
    for (int i=0; i<size; i++) {
      data[i] += a.data[i];
    }
  }
  void operator-=(const Vector &a) {
    assert(size == a.size);
    for (int i=0; i<size; i++) {
      data[i] -= a.data[i];
    }
  }
  void operator*=(double a) {
    for (int i=0; i<size; i++) {
      data[i] *= a;
    }
  }
private:
  int size;
  double *data;
};

inline Vector operator*(double a, const Vector &b) {
  const int sz = b.size;
  Vector out(sz);
  for (int i=0; i<sz; i++) {
    out.data[i] = a*b.data[i];
  }
  return out;
}
