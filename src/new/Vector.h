// -*- mode: C++; -*-

#pragma once

#include <cassert>
#include <string.h>
#include <math.h>

// A Vector is a reference-counted array of doubles.  You need to be
// careful, because a copy of a Vector (or the use of assignment,
// operator= when the array being assigned hasn't been initialized)
// will point to the same data as the original Vector, so if you want
// to create a new array that is a copy of an existing one, you should
// first construct the target array, and then do the assignment
// (i.e. operator=).  This is more than a little hokey, and I'm not
// sure whether it will lead to bugs.  My hope is that we won't often
// be *wanting* to copy an array.  After all, what's the point?

// One nice feature of Vector is that you can use the slice method to
// create a new vector pointing to a subset of the same data as an old
// one.  Once again, this allows shared data, hopefully enabling nice
// code while avoiding bloated memory use.

// Most arithmetic operators are defined on Vectors.  If you come
// across one that isn't defined, we could probably add its
// definition.

class Vector {
public:
  Vector() : size(0), offset(0), data(0), references_count(0) {}
  explicit Vector(int sz) : size(sz), offset(0), data(new double[size]), references_count(new int) {
    *references_count = 1;
  }
  Vector(const Vector &a) : size(a.size), offset(a.offset),
                            data(a.data), references_count(a.references_count) {
    *references_count += 1;
  }
  ~Vector() { free(); }
  void free() {
    if (references_count && *references_count) {
      *references_count -= 1;
      if (*references_count == 0) {
        delete[] data;
        delete references_count;
      }
      references_count = 0;
      data = 0;
      size = 0;
      offset = 0;
    }
  }
  void operator=(const Vector &a) {
    // I should note that the semantics of assignment are a little
    // tricky.
    if (!references_count) {
      size = a.size;
      offset = a.offset;
      data = a.data;
      references_count = a.references_count;
      *references_count += 1;
    } else {
      assert(size = a.size);
      memcpy(data+offset, a.data+a.offset, size*sizeof(double)); // faster than manual loop?
    }
  }

  // operator[] is a pair of "safe" operator for indexing a vector
  double operator[](int i) const {
    assert(i + offset >= 0);
    assert(i + offset < size);
    return data[i + offset];
  };
  double &operator[](int i) {
    assert(i + offset >= 0);
    assert(i + offset < size);
    return data[i + offset];
  };
  Vector slice(int start, int num) const {
    assert(start + num <= size);
    Vector out;
    out.size = num;
    out.data = data;
    out.offset = start;
    out.references_count = references_count;
    (*references_count) += 1;
    return out;
  }

  // Below are the arithmetic operators.  These should be made fast.
  Vector operator+(const Vector &a) const {
    // operator+ treats an empty vector as a zero vector;
    if (a.size == 0) return *this;
    if (size == 0) return a;
    assert(size == a.size);
    Vector out(size);
    const double *p1 = data + offset, *p2 = a.data + a.offset;
    for (int i=0; i<size; i++) {
      out.data[i] = p1[i] + p2[i];
    }
    return out;
  }
  Vector operator-(const Vector &a) const {
    // operator- treats an empty vector as a zero vector;
    if (a.size == 0) return *this;
    if (size == 0) return -a;
    assert(size == a.size);
    Vector out(size);
    const double *p1 = data + offset, *p2 = a.data + a.offset;
    for (int i=0; i<size; i++) {
      out.data[i] = p1[i] - p2[i];
    }
    return out;
  }
  Vector operator-() const {
    Vector out(size);
    const double *p1 = data + offset;
    for (int i=0; i<size; i++) {
      out.data[i] = -p1[i];
    }
    return out;
  }
  Vector operator*(double a) const {
    Vector out(size);
    const double *p2 = data + offset;
    for (int i=0; i<size; i++) {
      out.data[i] = a*p2[i];
    }
    return out;
  }
  friend Vector operator*(double a, const Vector &b);

  void operator+=(const Vector &a) {
    assert(size == a.size);
    double *p1 = data + offset;
    const double *p2 = a.data + a.offset;
    for (int i=0; i<size; i++) {
      p1[i] += p2[i];
    }
  }
  void operator-=(const Vector &a) {
    assert(size == a.size);
    double *p1 = data + offset;
    const double *p2 = a.data + a.offset;
    for (int i=0; i<size; i++) {
      p1[i] -= p2[i];
    }
  }
  void operator*=(double a) {
    double *p1 = data + offset;
    for (int i=0; i<size; i++) {
      p1[i] *= a;
    }
  }
  double sum() const {
    const double *p1 = data + offset;
    double out = 0;
    for (int i=0; i<size; i++) {
      out += p1[i];
    }
    return out;
  }
  double dot(const Vector &a) const {
    assert(a.size == size);
    double out = 0;
    const double *p1 = data + offset, *p2 = a.data + a.offset;
    for (int i=0; i<size; i++) {
      out += p1[i]*p2[i];
    }
    return out;
  }
  double norm() const {
    double out = 0;
    const double *p1 = data + offset;
    for (int i=0; i<size; i++) {
      out += p1[i]*p1[i];
    }
    return sqrt(out);
  }
  int get_size() const {
    return size;
  }
private:
  int size, offset;
  double *data;
  int *references_count; // counts how many objects refer to the data.
};

inline Vector operator*(double a, const Vector &b) {
  const int sz = b.size;
  Vector out(sz);
  const double *p2 = b.data + b.offset;
  for (int i=0; i<sz; i++) {
    out.data[i] = a*p2[i];
  }
  return out;
}
