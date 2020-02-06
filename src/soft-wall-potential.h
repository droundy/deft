#pragma once

// We assume z=0 to be overlap with the wall (and infinite potential)
static inline double soft_wall_potential(double z, double Vcutoff=1e300) {
  const double R_0 = 2*sigma/pow(2,5.0/6.0);
  const double rho = 1.0; // wall density

  if (z >= R_0) return 0;
  if ( z < 0 ) return Vcutoff;

  const double sig6 = uipow(sigma,6);
  const double sig12 = uipow(sigma,12);
  const double z3 = uipow(z,3);
  const double z9 = uipow(z,9);
  const double R3 = uipow(R_0,3);
  const double R9 = uipow(R_0,9);

  double potential = 2*M_PI*rho*epsilon*((z3-R3)/6
                                         + 2*sig12*(1/z9 - 1/R9)/45
                                         + (R_0 - z)*(R_0*R_0/2 + sig6/pow(R_0,4)
                                             - 2*sig12/5/pow(R_0,10))
                                         + sig6*(1/R3-1/z3)/3);
  if (potential < Vcutoff) return potential;
  return Vcutoff;
}
