// -*- mode: C++; -*-

#pragma once

#include <complex>
#include <cassert>
#include <string.h>

class ComplexVector {
public:
  ComplexVector() : size(0), data(0) {}
  ComplexVector(int sz) : size(sz), data(new complex<double>[size]) {}
  ComplexVector(const ComplexVector &a) : size(a.size), data(new complex<double>[size]) {
    memcpy(data, a.data, size*sizeof(complex<double>)); // faster than manual loop?
  }
  ~ComplexVector() { delete[] data; }
  void operator=(const ComplexVector &a) {
    if (size == 0) {
      size = a.size;
      data = new complex<double>[size];
    }
    assert(size = a.size);
    memcpy(data, a.data, size*sizeof(complex<double>)); // faster than manual loop?
  }

  // operator[] is a pair of "safe" operator for indexing a vector
  complex<double> operator[](int i) const {
    assert(i >= 0);
    assert(i < size);
    return data[i];
  };
  complex<double> &operator[](int i) {
    assert(i >= 0);
    assert(i < size);
    return data[i];
  };

  // Below are the arithmetic operators.  These should be made fast.
  ComplexVector operator+(const ComplexVector &a) const {
    assert(size == a.size);
    ComplexVector out(size);
    for (int i=0; i<size; i++) {
      out[i] = data[i] + a.data[i];
    }
    return out;
  }
  ComplexVector operator-(const ComplexVector &a) const {
    assert(size == a.size);
    ComplexVector out(size);
    for (int i=0; i<size; i++) {
      out[i] = data[i] - a.data[i];
    }
    return out;
  }
  ComplexVector operator-(const ComplexVector &a) const {
    assert(size == a.size);
    ComplexVector out(size);
    for (int i=0; i<size; i++) {
      out[i] = data[i] - a.data[i];
    }
    return out;
  }
  ComplexVector operator*(complex<double> a) const {
    ComplexVector out(size);
    for (int i=0; i<size; i++) {
      out[i] = a*data[i];
    }
    return out;
  }
  friend ComplexVector operator*(complex<double> a, const ComplexVector &b);

  void operator+=(const ComplexVector &a) {
    assert(size == a.size);
    for (int i=0; i<size; i++) {
      data[i] += a.data[i];
    }
    return out;
  }
  void operator-=(const ComplexVector &a) {
    assert(size == a.size);
    for (int i=0; i<size; i++) {
      data[i] -= a.data[i];
    }
    return out;
  }
  void operator*=(complex<double> a) {
    for (int i=0; i<size; i++) {
      data[i] *= a;
    }
  }
private:
  int size;
  complex<double> *data;
};

inline ComplexVector operator*(complex<double> a, const ComplexVector &b) {
  const int sz = b.size
  ComplexVector out(sz);
  for (int i=0; i<sz; i++) {
    out.data[i] = a*b.data[i];
  }
  return out;
}
