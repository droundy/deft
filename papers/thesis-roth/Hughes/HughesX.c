#include <stdio.h>
#include <math.h>

const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin
const double angstrom = 1.8897261; // An angstrom in atomic units

double hughes_X(double x) {
  const double R = 3.03420*angstrom/2;
  const double kT = 300*kB;
  const double epsilon_association = 1400.00*kB;
  const double epsilon_dispersion = 250.0*kB;
  const double lambda_dispersion = 1.78890;
  const double kappa_association = 0.0381876*(3.03420*angstrom*3.03420*angstrom*3.03420*angstrom);
  double output = 0;
  double boltz = -1.0*1.0 + exp(epsilon_association/kT);
  double n = x;
  double step = 4.188790204786391*R*R*R;
  double n3 = n*step;
  double deltak = 12.566370614359172*R*R;
  double n2 = deltak*n;
  double ghsyuwu = (R*n2*(5.555555555555555e-2*R*n2/(1.0 + -1.0*n3) + 0.5)/(1.0 + -1.0*n3) + 1.0)/(1.0 + -1.0*n3);
  double eta_d = 4.1887902047863905*R*R*R*n;
  double c1 = lambda_dispersion*(-1.50349*1.0 + 0.249434*lambda_dispersion) + 2.25855;
  double c2 = lambda_dispersion*(1.40049*1.0 + -0.827739*lambda_dispersion) + -0.66927;
  double c3 = lambda_dispersion*(-15.0427*1.0 + 5.30827*lambda_dispersion) + 10.1576;
  double eta_eff = eta_d*(eta_d*(c3*eta_d + c2) + c1);
  double ghs = (1.0 + -0.5*eta_eff)/((1.0 + -1.0*eta_eff)*(1.0 + -1.0*eta_eff)*(1.0 + -1.0*eta_eff));
  double deta_eff_by_deta_d = eta_d*(3.0*c3*eta_d + 2.0*c2) + c1;
  double dghs_by_deta_d = deta_eff_by_deta_d*((3.0*1.0 + -1.5*eta_eff)/(1.0 + -1.0*eta_eff) + -0.5)/((1.0 + -1.0*eta_eff)*(1.0 + -1.0*eta_eff)*(1.0 + -1.0*eta_eff));
  double da1_by_deta_d = epsilon_dispersion*(lambda_dispersion*lambda_dispersion*lambda_dispersion*(-4.0*dghs_by_deta_d*eta_d + -4.0*ghs) + 4.0*dghs_by_deta_d*eta_d + 4.0*ghs);
  double dc1_by_dlambda_dispersion = -1.50349*1.0 + 0.498868*lambda_dispersion;
  double dc2_by_dlambda_dispersion = 1.40049*1.0 + -1.655478*lambda_dispersion;
  double dc3_by_dlambda_dispersion = -15.0427*1.0 + 10.61654*lambda_dispersion;
  double deta_eff_by_dlambda_dispersion = eta_d*(eta_d*(dc3_by_dlambda_dispersion*eta_d + dc2_by_dlambda_dispersion) + dc1_by_dlambda_dispersion);
  double dghs_by_deta_eff = ((3.0*1.0 + -1.5*eta_eff)/(1.0 + -1.0*eta_eff) + -0.5)/((1.0 + -1.0*eta_eff)*(1.0 + -1.0*eta_eff)*(1.0 + -1.0*eta_eff));
  double dghs_by_dlambda_dispersion = deta_eff_by_dlambda_dispersion*dghs_by_deta_eff;
  printf("dghs_by_deta_eff is %g\n",dghs_by_deta_eff);
  double da1_by_dlambda_dispersion = epsilon_dispersion*eta_d*(lambda_dispersion*lambda_dispersion*(-4.0*dghs_by_dlambda_dispersion*lambda_dispersion + -12.0*ghs) + 4.0*dghs_by_dlambda_dispersion);
  double gSW = (-8.333333333333333e-2*da1_by_dlambda_dispersion*lambda_dispersion/eta_d + 0.25*da1_by_deta_d)/kT + ghsyuwu;
  printf("deta_eff_by_dlambda_dispersion is %g\n", deta_eff_by_dlambda_dispersion);
  printf("dghs_by_dlambda_dispersion %g\n", dghs_by_dlambda_dispersion);
  printf("da1 wrt eta is %g\n", da1_by_deta_d/kT);
  printf("da1 wrt lambda is %g\n", da1_by_dlambda_dispersion/kT);
  printf("gHS is %g\n", ghsyuwu);
  printf("gHS-effective is %g\n",ghs);
  printf("gSW is %g\n", gSW);
  double deltasaft = boltz*gSW*kappa_association;
  printf("Delta is %g\n", deltasaft);
  double n0 = 7.957747154594767e-2*deltak*n/(R*R);
  double X = (0.25*sqrt(8.0*deltasaft*n0 + 1.0) + -0.25)/(deltasaft*n0);
  output = X;
  
  return output;
}

void main() {
  const double n_cgs = 3.34e22; /*  atoms cm^-3 */
  const double n = 4.9388942e-3;
  printf("X is %g\n", hughes_X(n));
}
