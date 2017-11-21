// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include <stdio.h>

#pragma GCC diagnostic push
#if __GCC__ > 5
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Geometry>
#pragma GCC diagnostic pop

class Grid;

// The ifft is defined by:

// f(r) = 1/(2pi)^3 \int f(k) exp(-i k*r) d3r

Grid ifft(const GridDescription &gd, VectorXcd *rg);
Grid ifft(const GridDescription &gd, const VectorXcd &rg);

class ReciprocalGrid : public VectorXcd {
public:
  explicit ReciprocalGrid(const GridDescription &);
  ReciprocalGrid(const ReciprocalGrid &x);

  // We need to define this for our object to work
  typedef Eigen::VectorXcd Base;
  template<typename OtherDerived>
  ReciprocalGrid &operator=(const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }

  Grid ifft() const {
    return ::ifft(gd, *this);
  }

  void MultiplyBy(double f(Reciprocal));
  void MultiplyBy(complex f(Reciprocal));

  complex operator()(int x, int y, int z) const {
    if (z > gd.NzOver2) z = gd.Nz - z;
    return (*this)[x*gd.NyNzOver2 + y*gd.NzOver2 + z];
  }
  complex &operator()(int x, int y, int z) {
    if (z > gd.NzOver2) z = gd.Nz - z;
    return (*this)[x*gd.NyNzOver2 + y*gd.NzOver2 + z];
  }
  complex operator()(const Reciprocal &r) const {
    return (*this)(gd.Lat.toRelativeReciprocal(r));
  }
  complex operator()(const RelativeReciprocal &) const;
  GridDescription description() const { return gd; }
private:
  GridDescription gd;
};
