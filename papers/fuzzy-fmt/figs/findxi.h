#pragma once

static inline double find_alpha(double temp) {
  const double sigma=1;  //sigma must be 1 - changing it would invalidate other equations in the program!
  const double epsilon=1;
  return sigma*pow(2/(1+sqrt((temp*log(2))/epsilon)),1.0/6);
}

//Find Xi(T)---------------------------------------------------
static inline double Vwca(double r) {
    return 4.0*(pow(r,-12.0) - pow(r,-6.0)) +1;
}

static inline double f_wca(double r, double T){   //WCA mayer function
    return exp(-Vwca(r)/T) - 1;
}

static inline double B2_wca(double T) {
    long i=0;
    long num_points =10000;
    double r;
    double rmax_wca = pow(2.0,1.0/6);   //rmax_wca=1.122462048
    double dr=rmax_wca/num_points;
    double f_sum=0;
    for (; i<num_points; i++) {
      r=dr*(i+0.5);
      f_sum = f_sum + (-0.5)*(4*M_PI*r*r*dr*f_wca(r, T)); 
    }
    return f_sum;
}

static inline double B2_erf(double Xi, double T) {
  double alpha = find_alpha(T);
  return (M_PI/3)*((pow(alpha, 3) + 1.5*alpha*pow(Xi,2))*(1+erf(alpha/Xi)) + 1/pow(M_PI,0.5)*(pow(alpha,2)*Xi + pow(Xi,3))*exp(-pow((alpha/Xi),2)));
}


static inline double find_Xi(double T) {
  static double last_T = 0.0;
  static double last_Xi = 0.0;
  if (last_T == T) {
    return last_Xi;
  }
  double B2wca = B2_wca(T);
  double xi_lo = 0;
  double xi_hi = 1;
  double xi_mid;
  // printf("B2wca=%g\n", B2wca);
  do {
    xi_mid = 0.5*(xi_hi + xi_lo);
    if (B2_erf(xi_mid, T) > B2wca) {
      xi_hi = xi_mid;
    }  else  {
      xi_lo = xi_mid;
    }
  } while (xi_hi - xi_lo > 0.000000001);
  // printf("B2_erf mid=%g\n",B2_erf(xi_mid, T));
  last_T = T;
  last_Xi = xi_mid;
  return xi_mid;
}
//END Find Xi(T)---------------------------------------------------

static inline HomogeneousSFMTFluid sfmt_homogeneous(double n, double T) {
  HomogeneousSFMTFluid hf;
  hf.Xi() = find_Xi(T);
  hf.sigma() = 1;
  hf.epsilon() = 1;   //energy constant in the WCA fluid
  hf.kT() = T;
  hf.n() = n;
  hf.mu() = 0;
  return hf;
}



static inline SFMTFluidVeff sfmt_inhomogeneous(double T, double Lx, double Ly, double Lz, double dx) {
  SFMTFluidVeff hf(Lx, Ly, Lz, dx);
  hf.Xi() = find_Xi(T);
  hf.sigma() = 1;
  hf.epsilon() = 1;   //energy constant in the WCA fluid
  hf.kT() = T;
  hf.mu() = 0;
  return hf;
}
