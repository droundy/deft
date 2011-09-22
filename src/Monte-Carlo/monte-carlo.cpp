#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Monte-Carlo/monte-carlo.h"
#include "MersenneTwister.h"

int main(int argc, char *argv[]){
  if(argc!=5){
    printf("Incorrect number of args: radius, # of spheres, iterations, filename\n");
    return 1;
  }
  const double R = atof(argv[1]);
  const int Nspheres = atoi(argv[2]);
  const long Niter = atol(argv[3]);
  const char *filename = argv[4];
  run(R,Nspheres,Niter*Nspheres, filename);
}

double ran(){
  static MTRand my_mtrand;
  return my_mtrand.randExc(); // which is the range of [0,1)
}

Vector3d ran3(){
  double x, y, r2;
  do{
   x = 2 * ran() - 1;
   y = 2 * ran() - 1;
   r2 = x * x + y * y;
  } while(r2 >= 1 || r2 == 0);
  double fac = sqrt(-2*log(r2)/r2);
  Vector3d out(x*fac,y*fac,0);
  do{
   x = 2 * ran() - 1;
   y = 2 * ran() - 1;
   r2 = x * x + y * y;
  } while(r2 >= 1 || r2 == 0);
  fac = sqrt(-2*log(r2)/r2);
  out[2]=x*fac;
  return out;
}
