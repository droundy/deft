#pragma once

//This program finds alpha and xi for a given temperature T (really kT).
//Xi is the variance of r about the peak of the derivative of the mayer function with respect to the radius r
//Alpha is found by setting the second virial coefficient computed with V_wca equal 
//to the second virial coefficient computed with V_erf.
//It is referenced by new-melting.cpp

const double epsilon = 1.0;
const double sigma = 1.0;  //sigma must be 1, changing it would invalidate other equations

static const double Rmax = 2*sigma/pow(2.0,5.0/6.0); // this is the distance for which the potential is zero beyond it.

static inline double Vlj(double r) {
  const double power_6 = uipow(sigma/r, 6);
  return 4*epsilon*(sqr(power_6) - power_6);
}
static inline double Vwca(double r) {
  if (r > Rmax) return 0;
  return Vlj(r) + epsilon;
}

static inline double Vwca_prime(double r) {
  const double power_7 = uipow(sigma/r, 7);
  return 4*epsilon*(-12*sqr(power_7)*r+6*power_7);
}

static inline double f_wca(double r, double T) {
  return exp(-Vwca(r)/T) - 1;
}

static inline double fp_wca(double r, double T) {
  return exp(-Vwca(r)/T)*(-Vwca_prime(r)/T);
}

static inline double mean_fprime_wca_radius(double T) {
  const int N = 10000; // 1000000;
  double f_sum_top=0;
  double f_sum_bottom=0;
  const double dr = Rmax/N;
  for (double r=dr/2; r < Rmax; r += dr) {
  f_sum_top += r*dr*fp_wca(r,T);
  f_sum_bottom += dr*fp_wca(r,T);
  }
  return f_sum_top/f_sum_bottom;
}

static inline double mean_fprime_wca_radius_squared(double T) {
  const int N = 10000; // 1000000;
  double f_sum_top=0;
  double f_sum_bottom=0;
  const double dr = Rmax/N;
  for (double r=dr/2; r < Rmax; r += dr) {
  f_sum_top += uipow(r,2)*dr*fp_wca(r,T);
  f_sum_bottom += dr*fp_wca(r,T);
  }
  return f_sum_top/f_sum_bottom;
}
static inline double variance_fprime_wca_radius(double T) {
  return mean_fprime_wca_radius_squared(T) - uipow(mean_fprime_wca_radius(T),2);
}

static double cached_T = -1.0;
static double cached_alpha = -1.0;
static double cached_Xi = -1.0;

static inline double find_Xi(double T) {
  if (cached_T == T) {
    return cached_Xi;
  }
  return sqrt(2*variance_fprime_wca_radius(T));
}

static inline double B2_wca(double T) {   
  const int N = 10000; // 1000000;
  double f_sum=0;
  const double dr = Rmax/N;
  const double beta = 1.0/T;
  for (double r_cur=dr/2; r_cur < Rmax; r_cur += dr) {
  f_sum += -4*M_PI*r_cur*r_cur*dr*(exp(-beta*Vwca(r_cur)) - 1);
  }
  return f_sum/2;
}

static inline double B2_erf(double alpha, double Xi, double T) {
  return (M_PI/3)*((uipow(alpha, 3) + 1.5*alpha*uipow(Xi,2))*(1+erf(alpha/Xi))
                   + 1/sqrt(M_PI)*(uipow(alpha,2)*Xi + uipow(Xi,3))*exp(-uipow((alpha/Xi),2)));
}

//Find alpha(T) by matching second virial coefficeints B2 for WCA and Erf potentials: 
static inline double find_alpha(double T) {
  if (cached_T == T) {
    return cached_alpha;
  }
  double Xi = find_Xi(T);
  double B2wca = B2_wca(T);
  double alpha_lo = 0;
  double alpha_hi = 3*sigma;
  double alpha_mid;
  do {
    alpha_mid = 0.5*(alpha_hi + alpha_lo);
    if (B2_erf(alpha_mid, Xi, T) > B2wca) {
      alpha_hi = alpha_mid;
    }  else  {
      alpha_lo = alpha_mid;
    }
  } while (alpha_hi - alpha_lo > 0.000000001);
  printf("\nXi at T=%g is %g\n", T, find_Xi(T));
  printf("alpha at T=%g is %g\n", T, alpha_mid);
  printf("our B2 = %g\n", B2_erf(alpha_mid, Xi, T));
  printf("correct B2 = %g\n\n", B2_wca(T));
  cached_T = T;
  cached_Xi = Xi;
  cached_alpha = alpha_mid;
  return alpha_mid;
}

static inline HomogeneousSFMTFluid sfmt_homogeneous(double n, double T) {
  HomogeneousSFMTFluid hf;
  hf.Xi() = find_Xi(T);
  hf.sigma() = sigma;
  hf.epsilon() = epsilon;   //energy constant in the WCA fluid
  hf.kT() = T;
  hf.n() = n;
  hf.mu() = 0;
  return hf;
}

static inline SFMTFluidVeff sfmt_inhomogeneous(double T, double Lx, double Ly, double Lz, double dx) {
  SFMTFluidVeff hf(Lx, Ly, Lz, dx);
  hf.Xi() = find_Xi(T);
  hf.sigma() = sigma;
  hf.epsilon() = epsilon;   //energy constant in the WCA fluid
  hf.kT() = T;
  hf.mu() = 0;
  return hf;
}
