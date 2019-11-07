#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>

//Run this program from the directory it is listed in
//with command ./xi-fromB2

double sigma=1;
double epsilon=1;

double rmax_wca = pow(2.0,1.0/6);   //rmax_wca=1.122462048

double alpha(double T) {
    return pow(sqrt((2.0/(1+sqrt(T*log(2))))),1.0/3);
}

double rmax_erf(double T) { 
    return 10*alpha(T);
}

double Vwca(double r) {
    return 4.0*(pow(r,-12.0) - pow(r,-6.0)) +1;
}

double f_wca(double r, double T){   //WCA mayer function
    return exp(-Vwca(r)/T) - 1;
}

double B2_wca(double T) {   
    long i=0;
    long num_points =10000;
    double r;
    double dr=rmax_wca/num_points;
    double f_sum=0;
    for (; i<num_points; i++) {
      r=dr*i+ 0.0000000000001;
      f_sum = f_sum + (-0.5)*(4*M_PI*r*r*dr*f_wca(r, T)); 
    }
    return f_sum;
}


double f_erf(double r, double xi, double T){   //erf mayer function
    return (0.5)*(erf((r-alpha(T))/(xi/sqrt(2)))-1);
}

double B2_erf(double xi, double T) {   
    long i=0;
    long num_points =10000;
    double r;
    double dr=rmax_erf(T)/num_points;
    double f_sum=0;
    for (; i<num_points; i++) {
      r=dr*i+ 0.0000000000001;
      f_sum = f_sum + (-0.5)*(4*M_PI*r*r*dr*f_erf(r, xi, T)); 
    }
    return f_sum;
}


double find_Xi(double T) {      
    double B2wca = B2_wca(T);
    double xi_lo = 0;
    double xi_hi = 1;
    double xi_mid;
    printf("B2wca=%g\n", B2wca);
    do {
      xi_mid = 0.5*(xi_hi + xi_lo);
      if (B2_erf(xi_mid, T) > B2wca) {
        xi_hi = xi_mid;
      }  else  {
        xi_lo = xi_mid;
      } 
    } while (xi_hi - xi_lo > 0.000000001); 
    printf("B2_erf mid=%g\n",B2_erf(xi_mid, T));  
    printf("B2_erf hi =%g, xi_hi=%g\n",B2_erf(xi_hi, T), xi_hi);
    printf("B2_erf lo =%g, xi_lo=%g\n",B2_erf(xi_hi, T), xi_lo);
    printf("B2_erf=%g, xi=0.0457\n",B2_erf(0.0457, T));
    return xi_mid;  
}

int main() {
    long i=1;
    double T;
    double Xi_at_T=0;
    printf("rmax_wca=%g\n",rmax_wca); 
    for (; i<12; i++) {
      T=0.1*i; 
      //printf("rmax_erf=%g\n",rmax_erf(T)); 
      printf("T=%g\n", T);   
      Xi_at_T = find_Xi(T);
      printf("Xi=%g\n\n", Xi_at_T);
    }   

    return 0;
}
