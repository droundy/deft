// -*- mode: C++; -*-

#pragma once

#include "GridDescription.h"
#include <stdio.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Geometry>
#pragma GCC diagnostic pop

static const double default_eps_size = 1000.0;

class ReciprocalGrid;

ReciprocalGrid fft(const GridDescription &gd, const VectorXd &g);

class Grid : public VectorXd {
public:
  explicit Grid(const GridDescription &);
  Grid(const Grid &x);
  template <typename OtherDerived>
  explicit Grid(const GridDescription &gd, const Eigen::MatrixBase<OtherDerived> &x);

  // We need to define this for our object to work
  typedef Eigen::VectorXd Base;
  template<typename OtherDerived>
  Grid &operator=(const Eigen::MatrixBase <OtherDerived>& other) {
    this->Base::operator=(other);
    return *this;
  }
  ReciprocalGrid fft() const;
  double operator()(int x, int y, int z) const {
    return (*this)[x*gd.NyNz + y*gd.Nz + z];
  }
  double &operator()(int x, int y, int z) {
    return (*this)[x*gd.NyNz + y*gd.Nz + z];
  }
  double operator()(const Cartesian &r) const {
    return (*this)(gd.Lat.toRelative(r));
  }
  double operator()(const Relative &r) const;
  
  void Set(double f(Cartesian));
  void epsSlice(const char *fname, Cartesian xmax, Cartesian ymax,
                Cartesian corner, int resolution) const;
  void epsNativeSlice(const char *fname,
                      Cartesian xmax, Cartesian ymax, Cartesian corner) const;
  void epsNative1d(const char *fname, Cartesian xmin, Cartesian xmax, double yscale = 1, double xscale = 1, const char *comment = 0) const;
  void Dump1D(const char *fname, Cartesian xmin, Cartesian xmax) const;
  void epsRadial1d(const char *fname, double rmin = 0, double rmax = 0, double yscale = 1, double rscale = 1, const char *comment = 0) const;
  void ShellProjection(const VectorXd &R, VectorXd *output) const;
  double integrate() const {
    double val = 0;
    for (int i=0; i<gd.NxNyNz; i++) {
      val += (*this)[i];
    }
    return val*gd.dvolume;
  }
  GridDescription description() const { return gd; }
private:
  GridDescription gd;
};

template <typename OtherDerived> Grid::Grid(const GridDescription &gdin, const Eigen::MatrixBase<OtherDerived> &x)
  : VectorXd(gdin.NxNyNz), gd(gdin) {
  (*this) = x;
}
