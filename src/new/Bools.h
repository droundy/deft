// -*- mode: C++; -*-

#pragma once

#include <cassert>

class Bools {
public:
  Bools() : size(0), data(0) {}
  Bools(int sz) : size(sz), data(new bool[size]) {}
  Bools(const Bools &a) : size(a.size), data(new bool[size]) {
    memcpy(data, a.data, size*sizeof(bool)); // faster than manual loop?
  }
  void operator=(const Bools &a) {
    if (size == 0) {
      size = a.size;
      data = new bool[size];
    }
    assert(size = a.size);
    memcpy(data, a.data, size*sizeof(bool)); // faster than manual loop?
  }

  // operator[] is a pair of "safe" operator for indexing a vector
  bool operator[](int i) const {
    assert(i >= 0);
    assert(i < size);
    return data[i];
  };
  bool &operator[](int i) {
    assert(i >= 0);
    assert(i < size);
    return data[i];
  };

  // Below are the arithmetic operators.  These should be made fast.
  Bools operator&(const Bools &a) const {
    assert(size == a.size);
    Bools out(size);
    for (int i=0; i<size; i++) {
      out[i] = data[i] && a.data[i];
    }
    return out;
  }
  Bools operator|(const Bools &a) const {
    assert(size == a.size);
    Bools out(size);
    for (int i=0; i<size; i++) {
      out[i] = data[i] || a.data[i];
    }
    return out;
  }

  void operator&=(const Bools &a) {
    assert(size == a.size);
    for (int i=0; i<size; i++) {
      data[i] = data[i] && a.data[i];
    }
  }
  void operator|=(bool a) {
    for (int i=0; i<size; i++) {
      data[i] = data[i] || a;
    }
  }

  int count() const {
    int num = 0;
    for (int i=0; i<size; i++) {
      // The following uses a hokey attempt to ensure that we don't
      // ever end up with bools that have a value that is not either 1
      // or 0.
      num += data[i] & 1;
    }
    return num;
  }
private:
  int size;
  bool *data;
};
