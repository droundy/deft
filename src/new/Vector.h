// -*- mode: C++; -*-

#pragma once

#include <cassert>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <stdio.h>

#include "ComplexVector.h"

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
  Vector(double x, double y, double z) : size(3), offset(0), data(new double[3]), references_count(new int) {
    *references_count = 1;
    data[0] = x;
    data[1] = y;
    data[2] = z;
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
    if (!references_count && a.size) {
      size = a.size;
      offset = a.offset;
      data = a.data;
      references_count = a.references_count;
      *references_count += 1;
    } else if (!references_count) {
      size = 0;
    } else {
      assert(size == a.size);
      memcpy(data+offset, a.data+a.offset, size*sizeof(double)); // faster than manual loop?
    }
  }
  void operator=(double x) {
    for (int i=0; i<size; i++) {
      data[i + offset] = x;
    }
  }

  // operator[] is a pair of "safe" operator for indexing a vector
  double operator[](int i) const {
    assert(i + offset >= 0);
    assert(i < size);
    return data[i + offset];
  };
  double &operator[](int i) {
    assert(i + offset >= 0);
    assert(i < size);
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
  Vector operator+(double a) const {
    if (size == 0) return *this;
    Vector out(size);
    const double *p1 = data + offset;
    for (int i=0; i<size; i++) {
      out.data[i] = p1[i] + a;
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

  Vector operator/(double a) const {
    Vector out(size);
    const double *p2 = data + offset;
    for (int i=0; i<size; i++) {
      out.data[i] = p2[i]/a;
    }
    return out;
  }

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
  void operator/=(double a) {
    double *p1 = data + offset;
    for (int i=0; i<size; i++) {
      p1[i] /= a;
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
  double index3d(int Nx, int Ny, int Nz, int x, int y, int z) const {
    return (*this)[x*Ny*Nz + y*Nz + z];
    //return (*this)[x + y*Nx + z*Nx*Ny];
  }
  double &index3d(int Nx, int Ny, int Nz, int x, int y, int z) {
    return (*this)[x*Ny*Nz + y*Nz + z];
    //return (*this)[x + y*Nx + z*Nx*Ny];
  }
  void dumpSliceX(const char *fname, int Nx, int Ny, int Nz, int x) const {
    FILE *f = fopen(fname, "w");
    for (int y= Ny/2; y<Ny; y++) {
      for (int z=Nz/2; z<Nz; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      for (int z=0; z<=Nz/2; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      fprintf(f, "\n");
    }
    for (int y=0; y<=Ny/2; y++) {
      for (int z=Nz/2; z<Nz; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      for (int z=0; z<=Nz/2; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
  void dumpSliceY(const char *fname, int Nx, int Ny, int Nz, int y) const {
    FILE *f = fopen(fname, "w");
    for (int x= Nx/2; x<Nx; x++) {
      for (int z=Nz/2; z<Nz; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      for (int z=0; z<=Nz/2; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      fprintf(f, "\n");
    }
    for (int x=0; x<=Nx/2; x++) {
      for (int z=Nz/2; z<Nz; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      for (int z=0; z<=Nz/2; z++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
  void dumpSliceZ(const char *fname, int Nx, int Ny, int Nz, int z) const {
    FILE *f = fopen(fname, "w");
    for (int x= Nx/2; x<Nx; x++) {
      for (int y=Ny/2; y<Ny; y++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      for (int y=0; y<=Ny/2; y++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      fprintf(f, "\n");
    }
    for (int x=0; x<=Nx/2; x++) {
      for (int y=Ny/2; y<Ny; y++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      for (int y=0; y<=Ny/2; y++) {
        fprintf(f, "%g\t", index3d(Nx, Ny, Nz, x, y, z));
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
private:
  int size, offset;
  double *data;
  int *references_count; // counts how many objects refer to the data.
  friend Vector ifft(int Nx, int Ny, int Nz, double dV, ComplexVector f);
  friend ComplexVector fft(int Nx, int Ny, int Nz, double dV, Vector f);
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

inline ComplexVector fft(int Nx, int Ny, int Nz, double dV, Vector f) {
  assert(!(Nx&1)); // We want an even number of grid points in each direction.
  assert(!(Ny&1)); // We want an even number of grid points in each direction.
  assert(!(Nz&1)); // We want an even number of grid points in each direction.
  ComplexVector out(Nx*Ny*(int(Nz)/2 + 1));
  fftw_plan p = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, (double *)f.data+f.offset, (fftw_complex *)out.data, FFTW_WISDOM_ONLY);
  if (!p) {
    // It seems that fftw has not yet done enough measurement to make
    // a plan without modifying its input, so we need to take a moment
    // to time things to ensure a fast FFT.
    double *r = (double *)fftw_malloc((Nx*Ny*Nz+1)*sizeof(double)); // scratch space
    // First we make and destroy a plan with optimally aligned memory access...
    fftw_destroy_plan(fftw_plan_dft_r2c_3d(Nx, Ny, Nz, r, (fftw_complex *)out.data, FFTW_MEASURE));
    // Then we do the same with a poorly aligned array, so we're prepared for anything!
    fftw_destroy_plan(fftw_plan_dft_r2c_3d(Nx, Ny, Nz, r+1, (fftw_complex *)out.data, FFTW_MEASURE));
    fftw_free(r);
    // Now we will create the plan we actually use.
    p = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, (double *)f.data+f.offset, (fftw_complex *)out.data, FFTW_WISDOM_ONLY);
  }
  fftw_execute(p);
  fftw_destroy_plan(p);
  out *= dV;
  return out;
}

inline Vector ifft(int Nx, int Ny, int Nz, double dV, ComplexVector f) {
  assert(!(Nx&1)); // We want an even number of grid points in each direction.
  assert(!(Ny&1)); // We want an even number of grid points in each direction.
  assert(!(Nz&1)); // We want an even number of grid points in each direction.
  // Allocate a scratch array, since FFTW always overwrites its input
  // when performing a c2r transform.
  fftw_complex *c = (fftw_complex *)fftw_malloc(Nx*Ny*(int(Nz)/2+2)*sizeof(fftw_complex));
  memcpy(c, f.data+f.offset, 2*f.size*sizeof(double)); // faster than manual loop?
  Vector out(Nx*Ny*Nz); // create output vector
  fftw_plan p = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, c, (double *)out.data, FFTW_WISDOM_ONLY);
  if (!p) {
    // We need measurements!
    p = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, c, (double *)out.data, FFTW_MEASURE);
    // Now recopy data, which was trashed above
    memcpy(c, f.data+f.offset, 2*f.size*sizeof(double)); // faster than manual loop?
  }
  fftw_execute(p);
  fftw_destroy_plan(p);
  fftw_free(c);
  out *= 1.0/(Nx*Ny*Nz*dV);
  return out;
}
