// -*- mode: C++; -*-

#pragma once

#include <complex>
#include <cassert>
#include <string.h>
#include <math.h>


// A ComplexVector is a reference-counted array of std::complex<double>s.
// You need to be careful, because a copy of a ComplexVector (or the
// use of assignment, operator= when the array being assigned hasn't
// been initialized) will point to the same data as the original
// ComplexVector, so if you want to create a new array that is a copy
// of an existing one, you should first construct the target array,
// and then do the assignment (i.e. operator=).  This is more than a
// little hokey, and I'm not sure whether it will lead to bugs.  My
// hope is that we won't often be *wanting* to copy an array.  After
// all, what's the point?

// One nice feature of ComplexVector is that you can use the slice
// method to create a new vector pointing to a subset of the same data
// as an old one.  Once again, this allows shared data, hopefully
// enabling nice code while avoiding bloated memory use.

// Most arithmetic operators are defined on ComplexVectors.  If you come
// across one that isn't defined, we could probably add its
// definition.

class Vector;

class ComplexVector {
public:
  ComplexVector() : size(0), offset(0), data(0), references_count(0) {}
  explicit ComplexVector(long sz) : size(sz), offset(0), data(new std::complex<double>[size]),
                                   references_count(new int) {
    *references_count = 1;
  }
  ComplexVector(const ComplexVector &a) : size(a.size), offset(a.offset),
                            data(a.data), references_count(a.references_count) {
    *references_count += 1;
  }
  ~ComplexVector() { free(); }
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
  void operator=(const ComplexVector &a) {
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
      memcpy(data+offset, a.data+a.offset, size*sizeof(std::complex<double>)); // faster than manual loop?
    }
  }
  void operator=(std::complex<double> x) {
    for (long i=0; i<size; i++) {
      data[i] = x;
    }
  }

  // operator[] is a pair of "safe" operator for indexing a vector
  std::complex<double> operator[](long i) const {
    assert(i + offset >= 0);
    assert(i < size);
    return data[i + offset];
  };
  std::complex<double> &operator[](long i) {
    assert(i + offset >= 0);
    assert(i < size);
    return data[i + offset];
  };
  ComplexVector slice(long start, long num) const {
    assert(start + num <= size);
    ComplexVector out;
    out.size = num;
    out.data = data;
    out.offset = start;
    out.references_count = references_count;
    (*references_count) += 1;
    return out;
  }

  // Below are the arithmetic operators.  These should be made fast.
  ComplexVector operator+(const ComplexVector &a) const {
    // operator+ treats an empty vector as a zero vector;
    if (a.size == 0) return *this;
    if (size == 0) return a;
    assert(size == a.size);
    ComplexVector out(size);
    const std::complex<double> *p1 = data + offset, *p2 = a.data + a.offset;
    for (long i=0; i<size; i++) {
      out.data[i] = p1[i] + p2[i];
    }
    return out;
  }
  ComplexVector operator-(const ComplexVector &a) const {
    // operator- treats an empty vector as a zero vector;
    if (a.size == 0) return *this;
    if (size == 0) return -a;
    assert(size == a.size);
    ComplexVector out(size);
    const std::complex<double> *p1 = data + offset, *p2 = a.data + a.offset;
    for (long i=0; i<size; i++) {
      out.data[i] = p1[i] - p2[i];
    }
    return out;
  }
  ComplexVector operator-() const {
    ComplexVector out(size);
    const std::complex<double> *p1 = data + offset;
    for (long i=0; i<size; i++) {
      out.data[i] = -p1[i];
    }
    return out;
  }
  ComplexVector operator*(std::complex<double> a) const {
    ComplexVector out(size);
    const std::complex<double> *p2 = data + offset;
    for (long i=0; i<size; i++) {
      out.data[i] = a*p2[i];
    }
    return out;
  }
  friend ComplexVector operator*(std::complex<double> a, const ComplexVector &b);

  void operator+=(const ComplexVector &a) {
    assert(size == a.size);
    std::complex<double> *p1 = data + offset;
    const std::complex<double> *p2 = a.data + a.offset;
    for (long i=0; i<size; i++) {
      p1[i] += p2[i];
    }
  }
  void operator-=(const ComplexVector &a) {
    assert(size == a.size);
    std::complex<double> *p1 = data + offset;
    const std::complex<double> *p2 = a.data + a.offset;
    for (long i=0; i<size; i++) {
      p1[i] -= p2[i];
    }
  }
  void operator*=(std::complex<double> a) {
    std::complex<double> *p1 = data + offset;
    for (long i=0; i<size; i++) {
      p1[i] *= a;
    }
  }
  std::complex<double> sum() const {
    const std::complex<double> *p1 = data + offset;
    std::complex<double> out = 0;
    for (long i=0; i<size; i++) {
      out += p1[i];
    }
    return out;
  }
  std::complex<double> dot(const ComplexVector &a) const {
    assert(a.size == size);
    std::complex<double> out = 0;
    const std::complex<double> *p1 = data + offset, *p2 = a.data + a.offset;
    for (long i=0; i<size; i++) {
      out += p1[i]*p2[i];
    }
    return out;
  }
  std::complex<double> norm() const {
    std::complex<double> out = 0;
    const std::complex<double> *p1 = data + offset;
    for (long i=0; i<size; i++) {
      out += p1[i]*p1[i];
    }
    return sqrt(out);
  }
  long get_size() const {
    return size;
  }
  std::complex<double> index3d(long Nx, long Ny, long Nz, long x, long y, long z) const {
    return (*this)[x*Ny*(Nz/2+1) + y*(Nz/2+1) + z];
  }
  std::complex<double> &index3d(long Nx, long Ny, long Nz, long x, long y, long z) {
    return (*this)[x*Ny*(Nz/2+1) + y*(Nz/2+1) + z];
  }
  void dumpSliceXR(const char *fname, long Nx, long Ny, long Nz, long x) const {
    FILE *f = fopen(fname, "w");
    if (!f) {
      printf("unable to create file %s\n", fname);
      return;
    }
    for (long y= Ny/2; y<Ny; y++) {
      for (long z=0; z<=Nz/2; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z).real());
      }
      fprintf(f, "\n");
    }
    for (long y=0; y<=Ny/2; y++) {
      for (long z=0; z<=Nz/2; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z).real());
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
  void dumpSliceXI(const char *fname, long Nx, long Ny, long Nz, long x) const {
    FILE *f = fopen(fname, "w");
    if (!f) {
      printf("unable to create file %s\n", fname);
      return;
    }
    for (long y= Ny/2; y<Ny; y++) {
      for (long z=0; z<=Nz/2; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z).imag());
      }
      fprintf(f, "\n");
    }
    for (long y=0; y<=Ny/2; y++) {
      for (long z=0; z<=Nz/2; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z).imag());
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }

private:
  long size, offset;
  std::complex<double> *data;
  int *references_count; // counts how many objects refer to the data.
  friend Vector ifft(long Nx, long Ny, long Nz, double dV, ComplexVector f);
  friend ComplexVector fft(long Nx, long Ny, long Nz, double dV, Vector f);
};

inline ComplexVector operator*(std::complex<double> a, const ComplexVector &b) {
  const long sz = b.size;
  ComplexVector out(sz);
  const std::complex<double> *p2 = b.data + b.offset;
  for (long i=0; i<sz; i++) {
    out.data[i] = a*p2[i];
  }
  return out;
}
