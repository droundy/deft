#pragma once

static inline double find_alpha(double temp) {
  assert(sigma == 1.0); //sigma must be 1 - changing it would invalidate other equations in the program?
  return sigma*pow(2/(1+sqrt((temp*log(2))/epsilon)),1.0/6);
}


//Old Xi derived from derivatives (referenced here if want to revert to old Xi for comparison)--------DEBUG
// static inline double find_Xi(double kT) {
//   double xi = find_alpha(kT)/(6*pow(M_PI,0.5))/(log(2) + pow(log(2)*epsilon/kT,0.5));  //check syntex!
//   printf("our alpha at T=%g is %g\n", kT, find_alpha(kT));
//   printf("our Xi at T=%g is %g\n", kT, xi);
//   return xi;
// }
//END old Xi(T)---------------------------------------------------------------------------------------DEBUG



//Find Xi(T)---------------------------------------------------
static inline double Vwca(double r) {
    return 4.0*epsilon*(pow(r/sigma,-12.0) - pow(r/sigma,-6.0)) + epsilon;
}

static inline double f_wca(double r, double T){   //WCA mayer function
    return exp(-Vwca(r)/T) - 1;
}

static inline double B2_wca(double T) {
    long num_points =10000;
    double Rmax = sigma*pow(2.0,1.0/6);   //rmax_wca=1.122462048
    double dr = Rmax/num_points;
    double f_sum=0;
    for (double r=dr/2; r < Rmax; r += dr) {
      f_sum = f_sum + (-0.5)*(4*M_PI*r*r*dr*f_wca(r, T));
    }
    return f_sum;
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
  //putStrLn("B2wca=%g\n", B2wca);   //DEBUG
  do {
    xi_mid = 0.5*(xi_hi + xi_lo);
    if (B2_erf(xi_mid, T) > B2wca) {
      xi_hi = xi_mid;
    }  else  {
      xi_lo = xi_mid;
    }
  } while (xi_hi - xi_lo > 0.000000001);
  //putStrLn("B2_erf mid=%g\n",B2_erf(xi_mid, T));   //DEBUG
  last_T = T;
  last_Xi = xi_mid;
  printf("our alpha at T=%g is %g\n", T, find_alpha(T));
  printf("our Xi at T=%g is %g\n", T, xi_mid);
  printf("our B2 = %g\n", B2_erf(xi_mid, T));
  printf("correct B2 = %g\n", B2_wca(T));
  return xi_mid;
}
//END Find Xi(T)---------------------------------------------------



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
