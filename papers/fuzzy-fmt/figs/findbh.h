#pragma once

const double epsilon = 1.0;
const double sigma = 1.0;  //sigma must be 1, changing it would invalidate other equations

//When edited new-bh-radial-wca.cpp and this file to be simliar to new-radial-wca.cpp and findxi.h, R and sigma were CHANGED!
//WAS: R = 2*radius, sigma=R*pow(2,-1.0/6.0)  NOW: R = 2*sigma/pow(2.0,5.0/6.0), sigma=1

static const double Rmax = 2*sigma/pow(2.0,5.0/6.0); // this is the distance for which the potential is zero beyond it.
const double radius = Rmax/2.0;

static inline double Vwca(double r) {
  const double power_6 = uipow(sigma/r, 6);
  if (r > Rmax) return 0;
  return 4*epsilon*(sqr(power_6) - power_6) + epsilon;
}

static inline double R_BH(const double T) {          //note: T means kT
  const int N = 10000; // 1000000; ask David which to use 10000 or 1000000 (original)
  double bh_diameter = 0;
  const double dr = Rmax/N;
  const double beta = 1.0/T;
  for (double r_cur=dr/2; r_cur < Rmax; r_cur += dr) {
    bh_diameter += (1 - exp(-beta*Vwca(r_cur)))*dr;
  }
  printf("barker_henderson diameter at T=%g is %g\n", T, bh_diameter);
  printf("barker_henderson B2 = %g\n", 2*M_PI/3*uipow(bh_diameter,3));
  return bh_diameter/2;
}

static inline HomogeneousWhiteBearFluid bh_homogeneous(double n, double T) {
  HomogeneousWhiteBearFluid hf;

  double rad_bh = R_BH(T);
  printf("rad_bh is %g, sigma is %g and packing fraction is %g\n", rad_bh, sigma, n*uipow(rad_bh, 3)*4*M_PI/3);

  hf.R() = rad_bh;
  hf.kT() = T;
  hf.n() = n;   //ASK! should this be replaced with line below??
  //hf.n() = reduced_density*pow(2,-5.0/2.0);   //OLD line in new-bh-radial-wca.cpp -ASK 
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
