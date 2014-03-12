//From ~/deft/src/WaterSaftFast.cpp

#include <stdio.h>
#include <math.h>

const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin
const double angstrom = 1.8897261; // An angstrom in atomic units
const double Ha = 4.395e-18; // Hartree in Joules
const double atm = 101325; // atm in Pascal
const double bohr = 5.29177208e-9; // bohr radius in cm
const double m = 18.01528; // molar mass of water, g mol^-1
const double N_A = 6.02213e23; // mol^-1

double hughes_Fsaft(double x)
{
  const double R = 3.03420*angstrom/2;
  const double kT = 300*kB;
  const double epsilon_association = 1400.00*kB;
  const double epsilon_dispersion = 250.0*kB;
  const double lambda_dispersion = 1.78890;
  const double kappa_association = 0.0381876*(3.03420*angstrom*3.03420*angstrom*3.03420*angstrom);  
  double output = 0;
  double deltak = 12.566370614359172*R*R;
  double n = x;
  double n0 = 7.957747154594767e-2*deltak*n/(R*R);
  double boltz = -1.0*1.0 + exp(epsilon_association/kT);
  double step = 4.188790204786391*R*R*R;
  double n3 = n*step;
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
  double da1_by_deta_d = epsilon_dispersion*(dghs_by_deta_d*eta_d + ghs)*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
  double dc1_by_dlambda_dispersion = -1.50349*1.0 + 0.498868*lambda_dispersion;
  double dc2_by_dlambda_dispersion = 1.40049*1.0 + -1.655478*lambda_dispersion;
  double dc3_by_dlambda_dispersion = -15.0427*1.0 + 10.61654*lambda_dispersion;
  double deta_eff_by_dlambda_dispersion = eta_d*(eta_d*(dc3_by_dlambda_dispersion*eta_d + dc2_by_dlambda_dispersion) + dc1_by_dlambda_dispersion);
  double dghs_by_dlambda_dispersion = deta_eff_by_dlambda_dispersion*((3.0*1.0 + -1.5*eta_eff)/(1.0 + -1.0*eta_eff) + -0.5)/((1.0 + -1.0*eta_eff)*(1.0 + -1.0*eta_eff)*(1.0 + -1.0*eta_eff));
  double da1_by_dlambda_dispersion = epsilon_dispersion*eta_d*(-12.0*ghs*lambda_dispersion*lambda_dispersion + dghs_by_dlambda_dispersion*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0));
  double gSW = (-8.333333333333333e-2*da1_by_dlambda_dispersion*lambda_dispersion/eta_d + 0.25*da1_by_deta_d)/kT + ghsyuwu;
  double deltasaft = boltz*gSW*kappa_association;
  double X = (0.25*sqrt(8.0*deltasaft*n0 + 1.0) + -0.25)/(deltasaft*n0);
  double Fassoc = kT*n0*(2.0*1.0 + 4.0*log(X) + -2.0*X);
  printf("Fassoc is %g\n", Fassoc);
  double a1 = epsilon_dispersion*eta_d*ghs*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
  double a1integrated = a1*n;
  double KHS = (eta_d*(eta_d*(eta_d*(-4.0*1.0 + eta_d) + 6.0) + -4.0) + 1.0)/(eta_d*(4.0*1.0 + 4.0*eta_d) + 1.0);
  double a2 = 0.5*KHS*da1_by_deta_d*epsilon_dispersion*eta_d;
  double a2integrated = a2*n/kT;
  double Fdisp = a2integrated + a1integrated;
  printf("Fdisp is %g\n", Fdisp);
  double Fideal = kT*n*(-1.0*1.0 + log(2.646476976618268e-6*n/(sqrt(kT)*kT)));
  printf("Fideal is %g\n", Fideal);
  double phi1 = -1.0*n0*log(1.0 + -1.0*n3);
  double kTphi1 = kT*phi1;
  double n1 = 7.957747154594767e-2*deltak*n/R;
  double phi2 = n1*n2/(1.0 + -1.0*n3);
  double kTphi2 = kT*phi2;
  double phi3 = n2*n2*n2*((log(1.0 + -1.0*n3)*(8.841941282883075e-3*1.0/n3 + -1.768388256576615e-2) + 8.841941282883075e-3)/n3 + 8.841941282883075e-3*log(1.0 + -1.0*n3))/((1.0 + -1.0*n3)*(1.0 + -1.0*n3));
  double kTphi3 = kT*phi3;
  double whitebear = kTphi3 + kTphi2 + kTphi1;
  printf("Fhs is %g\n", whitebear);
  const double mu=0; //My Python code does not invovle mu
  double FSAFT = mu*n + whitebear + Fideal + Fdisp + Fassoc;
  output = FSAFT;

  return output;

}

double dF_dn(double n)
{
  const double dn = 1e-10;
  printf("dF_dn is %g\n",(hughes_Fsaft(n+dn) - hughes_Fsaft(n))/dn);

  return (hughes_Fsaft(n+dn) - hughes_Fsaft(n))/dn;
}

double pressure(double n)
{
  return n*dF_dn(n) - hughes_Fsaft(n);
}

void main()
{
  const double n = 4.9388942e-3; 
  // const double n = 1e-19;
  const double conv_p = Ha*1e6/pow(bohr,3.0)/atm;
  const double conv_n = m/pow(bohr,3.0)/N_A;
  printf("n in g/mL is %g\n",n*conv_n);
  printf("Pressure in atm is %g\n",pressure(n)*conv_p);
  printf("F is %g\n", hughes_Fsaft(n));
  printf("F with small change is %g\n", hughes_Fsaft(n+1e-10));
}
