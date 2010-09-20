// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include <stdio.h>
#include <Eigen/Geometry>

class Grid;

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
private:
  GridDescription gd;
};
