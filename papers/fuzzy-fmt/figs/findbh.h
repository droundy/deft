#pragma once

static const double Rmax = 2*sigma/pow(2.0,5.0/6.0); // this is the distance for which the potential is zero beyond it.

static inline double Vwca(double r) {
  const double power_6 = uipow(sigma/r, 6);
  if (r > Rmax) return 0;
  return 4*epsilon*(sqr(power_6) - power_6) + epsilon;
}

static inline double R_BH(const double kT) {
  const int N = 10000; // 1000000;
  double bh_diameter = 0;
  const double dr = Rmax/N;
  const double beta = 1.0/kT;
  for (double r_cur=dr/2; r_cur < Rmax; r_cur += dr) {
    bh_diameter += (1 - exp(-beta*Vwca(r_cur)))*dr;
  }
  printf("barker_henderson diameter at T=%g is %g\n", kT, bh_diameter);
  printf("barker_henderson B2 = %g\n", 2*M_PI/3*uipow(bh_diameter,3));
  return bh_diameter/2;
}

static inline HomogeneousWhiteBearFluid bh_homogeneous(double n, double T) {
  HomogeneousWhiteBearFluid hf;

  double rad_bh = R_BH(T);
  printf("rad_bh is %g, sigma is %g and packing fraction is %g\n", rad_bh, sigma, n*uipow(rad_bh, 3)*4*M_PI/3);

  hf.R() = rad_bh;
  hf.kT() = T;
  hf.n() = n;
  hf.mu() = 0;
  return hf;
}

static inline WhiteBearFluidVeff bh_inhomogeneous(double T, double Lx, double Ly, double Lz, double dx) {
  WhiteBearFluidVeff f(Lx, Ly, Lz, dx);
  f.R() = R_BH(T);
  f.kT() = T;
  f.mu() = 0;
  return f;
}
