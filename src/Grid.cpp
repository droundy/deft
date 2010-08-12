#include "Grid.h"
#include "ReciprocalGrid.h"
#include <fftw3.h>

double Grid::operator()(const Relative &r) const {
  double rx = r(0)*gd.Nx, ry = r(1)*gd.Ny, rz = r(2)*gd.Nz;
  int ix = floor(rx), iy = floor(ry), iz = floor(rz);
  double wx = rx-ix, wy = ry-iy, wz = rz-iz;
  while (ix < 0) ix += gd.Nx;
  while (iy < 0) iy += gd.Ny;
  while (iz < 0) iz += gd.Nz;
  while (ix >= gd.Nx) ix -= gd.Nx;
  while (iy >= gd.Ny) iy -= gd.Ny;
  while (iz >= gd.Nz) iz -= gd.Nz;
  int ixp1 = (ix+1)%gd.Nx, iyp1 = (iy+1)%gd.Ny, izp1 = (iz+1)%gd.Nz;
  while (ixp1 < 0) ixp1 += gd.Nx;
  while (iyp1 < 0) iyp1 += gd.Ny;
  while (izp1 < 0) izp1 += gd.Nz;
  while (ixp1 >= gd.Nx) ixp1 -= gd.Nx;
  while (iyp1 >= gd.Ny) iyp1 -= gd.Ny;
  while (izp1 >= gd.Nz) izp1 -= gd.Nz;
  assert(wx>=0);
  assert(wy>=0);
  assert(wz>=0);
  return (1-wx)*(1-wy)*(1-wz)*(*this)(ix,iy,iz)
    + wx*(1-wy)*(1-wz)*(*this)(ixp1,iy,iz)
    + (1-wx)*wy*(1-wz)*(*this)(ix,iyp1,iz)
    + (1-wx)*(1-wy)*wz*(*this)(ix,iy,izp1)
    + wx*(1-wy)*wz*(*this)(ixp1,iy,izp1)
    + (1-wx)*wy*wz*(*this)(ix,iyp1,izp1)
    + wx*wy*(1-wz)*(*this)(ixp1,iyp1,iz)
    + wx*wy*wz*(*this)(ixp1,iyp1,izp1);
}

Grid::Grid(const GridDescription &gdin)
  : VectorXd(gdin.NxNyNz), gd(gdin),
    r2_op(gd, cartSqr),
    x_op(gd, xfunc),
    y_op(gd, yfunc),
    z_op(gd, zfunc) {
}

Grid::Grid(const Grid &x) : VectorXd(x), gd(x.gd),
                            r2_op(gd, cartSqr),
                            x_op(gd, xfunc), y_op(gd, yfunc), z_op(gd, zfunc) {
}

void Grid::Set(double f(Cartesian)) {
  for (int x=0; x<gd.Nx; x++) {
    for (int y=0; y<gd.Ny; y++) { 
      for (int z=0; z<gd.Nz; z++) {
        Cartesian h(gd.Lat.wignerSeitz(gd.Lat.toCartesian(Relative(x*gd.dx,y*gd.dy,z*gd.dz))));
        (*this)(x,y,z) = f(h);
      }
    }
  }
}

void Grid::epsSlice(const char *fname,
                    Cartesian xmax, Cartesian ymax, Cartesian corner,
                    int resolution) const {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
  } else {
    fprintf(out, "%%!PS-Adobe-3.0 EPSF\n");
    const double size = xmax.norm() + ymax.norm();
    Cartesian ddx(xmax/resolution);
    Cartesian ddy(ymax/resolution);
    Cartesian xhat(ddx), yhat(ddy);
    xhat /= xhat.norm();
    yhat -= xhat*yhat.dot(xhat);
    yhat /= yhat.norm();
    fprintf(out, "%%%%BoundingBox: 0 0 %lg %lg\n",
            xmax.norm()*default_eps_size/size,
            ymax.norm()*default_eps_size/size);
    fprintf(out, "gsave\n");
    fprintf(out, "%lg %lg scale\n", default_eps_size/size,
            default_eps_size/size);
    //fprintf(out, "%lg %lg translate\n", -xmin, -ymin);
    fprintf(out, "/Times-Roman findfont 20 scalefont setfont\n");
    fprintf(out, "newpath 140 280 moveto (%s) show\n", fname);
    double max = this->maxCoeff();
    if (-this->minCoeff() > this->maxCoeff()) {
      max = -this->minCoeff();
    }
    fprintf(out, "/max %lg def\n", max);
    fprintf(out, "1 setlinecap\n");
    fprintf(out, "/P {\n\
    max div\n\
    dup 0 lt {\n\
        1 add\n\
        dup 1\n\
    }{\n\
        neg 1 add\n\
        dup 1 3 1 roll\n\
    } ifelse\n\
    setrgbcolor\n\
    newpath\n\
    moveto\n");
    fprintf(out, "    %g %g rmoveto\n",
            (-corner + (ddx+ddy)*0.5).dot(xhat),
            (-corner + (ddx+ddy)*0.5).dot(yhat));
    fprintf(out, "    %g %g rlineto\n", -ddy.dot(xhat), -ddy.dot(yhat));
    fprintf(out, "    %g %g rlineto\n", -ddx.dot(xhat), -ddx.dot(yhat));
    fprintf(out, "    %g %g rlineto\n", ddy.dot(xhat), ddy.dot(yhat));
    fprintf(out, "    %g %g rlineto\n", ddx.dot(xhat), ddx.dot(yhat));
    fprintf(out, "    gsave\n\
    fill\n                     \
    grestore\n");
    fprintf(out, "    %g setlinewidth\n", 0.1*ddx.norm());
    fprintf(out, "    stroke\n\
} def\n");

    // We now just need to output the actual data!
    for (int x=0; x<=resolution; x++) {
      for (int y=0; y<=resolution; y++) {
        Cartesian here(corner + x*ddx + y*ddy);
        const double fhere = (*this)(here);
        //fprintf(out, "%% %g\t%g\t%g\n", here(0), here(1), here(2));
        fprintf(out, "%g\t%g\t%g\tP\n",
                x*ddx.norm() + corner.dot(xhat),
                y*ddy.norm() + corner.dot(yhat), fhere);
      }
    }
    
    // And now we can output the trailer!
    fprintf(out, "grestore\n");
    fprintf(out, "showpage\n");
    fprintf(out, "%%%%Trailer\n");
    fprintf(out, "%%%%EOF\n");
    fclose(out);
  }
}

void Grid::epsNativeSlice(const char *fname,
                          Cartesian xmax, Cartesian ymax,
                          Cartesian corner) const {
  FILE *out = fopen(fname, "w");
  if (!out) {
    fprintf(stderr, "Unable to create file %s!\n", fname);
    // don't just abort?
  } else {
    // put corner on the grid!
    corner = gd.fineLat.round(corner);
    fprintf(out, "%%!PS-Adobe-3.0 EPSF\n");
    const double size = xmax.norm() + ymax.norm();
    Lattice weelat(gd.fineLat);
    weelat.reorientBasis(Cartesian(ymax.cross(xmax)));
    Cartesian ddx(weelat.a1());
    Cartesian ddy(weelat.a2());
    Cartesian xhat(xmax), yhat(ymax);
    xhat /= xhat.norm();
    yhat -= xhat*yhat.dot(xhat);
    yhat /= yhat.norm();
    fprintf(out, "%%%%BoundingBox: 0 0 %lg %lg\n",
            xmax.norm()*default_eps_size/size,
            ymax.norm()*default_eps_size/size);
    fprintf(out, "%% ddx is %g %g %g\n", ddx(0), ddx(1), ddx(2));
    fprintf(out, "gsave\n");
    fprintf(out, "%lg %lg scale\n", default_eps_size/size,
            default_eps_size/size);
    //fprintf(out, "%lg %lg translate\n", -xmin, -ymin);
    fprintf(out, "/Times-Roman findfont 20 scalefont setfont\n");
    fprintf(out, "newpath 140 280 moveto (%s) show\n", fname);
    double max = this->maxCoeff();
    if (-this->minCoeff() > this->maxCoeff()) {
      max = -this->minCoeff();
    }
    fprintf(out, "/max %lg def\n", max);
    fprintf(out, "1 setlinecap\n");
    fprintf(out, "/P {\n\
    max div\n\
    dup 0 lt {\n\
        1 add\n\
        dup 1\n\
    }{\n\
        neg 1 add\n\
        dup 1 3 1 roll\n\
    } ifelse\n\
    setrgbcolor\n\
    newpath\n\
    moveto\n");
    fprintf(out, "    %g %g rmoveto\n",
            (-corner + (ddx+ddy)*0.5).dot(xhat),
            (-corner + (ddx+ddy)*0.5).dot(yhat));
    fprintf(out, "    %g %g rlineto\n", -ddy.dot(xhat), -ddy.dot(yhat));
    fprintf(out, "    %g %g rlineto\n", -ddx.dot(xhat), -ddx.dot(yhat));
    fprintf(out, "    %g %g rlineto\n", ddy.dot(xhat), ddy.dot(yhat));
    fprintf(out, "    %g %g rlineto\n", ddx.dot(xhat), ddx.dot(yhat));
    fprintf(out, "    gsave\n\
    fill\n                     \
    grestore\n");
    fprintf(out, "    %g setlinewidth\n", 0.1*ddx.norm());
    fprintf(out, "    stroke\n\
} def\n");

    // We now just need to output the actual data!
    const double small = ddx.norm() + ddy.norm();
    const int nmax = 4*(ymax.norm() + xmax.norm())/small;
    for (int x=-nmax; x<=nmax; x++) {
      for (int y=-nmax; y<=nmax; y++) {
        Cartesian dr(x*ddx + y*ddy);
        Cartesian here(corner + x*ddx + y*ddy);
        const double fhere = (*this)(here);
        if (dr.dot(xhat) > -small &&
            dr.dot(xhat) < xmax.dot(xhat) + small &&
            dr.dot(yhat) > -small &&
            dr.dot(yhat) < ymax.dot(yhat) + small) {
          //fprintf(out, "%% %g\t%g\t%g\n", here(0), here(1), here(2));
          fprintf(out, "%g\t%g\t%g\tP\n",
                  here.dot(xhat), here.dot(yhat), fhere);
        }
      }
    }
    
    // And now we can output the trailer!
    fprintf(out, "grestore\n");
    fprintf(out, "showpage\n");
    fprintf(out, "%%%%Trailer\n");
    fprintf(out, "%%%%EOF\n");
    fclose(out);
  }
}

ReciprocalGrid Grid::fft() const {
  ReciprocalGrid out(gd);
  const double *mydata = data();
  fftw_execute(fftw_plan_dft_r2c_3d(gd.Nx, gd.Ny, gd.Nz, (double *)mydata, (fftw_complex *)out.data(), FFTW_ESTIMATE));
  out /= gd.NxNyNz;
  return out;
}
