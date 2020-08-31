#pragma once

const double epsilon = 1.0;
const double sigma = 1.0;  //sigma must be 1, changing it would invalidate other equations

static inline double find_alpha(double T) {                    //note: T means kT
  return sigma*pow(2/(1+sqrt((T*log(2))/epsilon)),1.0/6);
}

static const double Rmax = 2*sigma/pow(2.0,5.0/6.0); // this is the distance for which the potential is zero beyond it.

static inline double Vlj(double r) {
  const double power_6 = uipow(sigma/r, 6);
  return 4*epsilon*(sqr(power_6) - power_6);
}
static inline double Vwca(double r) {
  if (r > Rmax) return 0;
  return Vlj(r) + epsilon;
}

//Find Xi(T) by matching second virial coefficeints B2 for WCA and Erf potentials: 

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

static inline double B2_erf(double Xi, double T) {
  double alpha = find_alpha(T);
  return (M_PI/3)*((uipow(alpha, 3) + 1.5*alpha*uipow(Xi,2))*(1+erf(alpha/Xi))
                   + 1/sqrt(M_PI)*(uipow(alpha,2)*Xi + uipow(Xi,3))*exp(-uipow((alpha/Xi),2)));
}

static inline double find_Xi(double T) {
  static double last_T = 0.0;
  static double last_Xi = 0.0;
  if (last_T == T) {
    return last_Xi;
  }
  double B2wca = B2_wca(T);
  double xi_lo = 0;
  double xi_hi = sigma;
  double xi_mid;
  do {
    xi_mid = 0.5*(xi_hi + xi_lo);
    if (B2_erf(xi_mid, T) > B2wca) {
      xi_hi = xi_mid;
    }  else  {
      xi_lo = xi_mid;
    }
  } while (xi_hi - xi_lo > 0.000000001);
  last_T = T;
  last_Xi = xi_mid;
  printf("our alpha at T=%g is %g\n", T, find_alpha(T));
  printf("our Xi at T=%g is %g\n", T, xi_mid);
  printf("our B2 = %g\n", B2_erf(xi_mid, T));
  printf("correct B2 = %g\n", B2_wca(T));
  return xi_mid;
}

//Old Xi derived from derivatives (referenced here if want to revert to old Xi for comparison)
 //static inline double find_Xi(double T) {
   //double xi = find_alpha(T)/(6*uipow(M_PI,0.5))/(log(2) + uipow(log(2)*epsilon/T,0.5));
   //printf("our alpha at T=%g is %g\n", T, find_alpha(T));
   //printf("our Xi at T=%g is %g\n", T, xi);
   //return xi;
 //}

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
