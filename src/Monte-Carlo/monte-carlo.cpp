#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Monte-Carlo/monte-carlo.h"


int main(){
  printf("Hello World\n");
  const int N = 14;
  Vector3d *spheres = new Vector3d[N];
  for(int i = 0; i < N; i++){
    //  for(int j = 0; j < 3; j++){
    //  spheres[i][j] = ran();
    // }
    spheres[i] = ran3();
  }
  for(int i = 0; i<N; i++){
    printf("%g  %g  %g\n", spheres[i][0], spheres[i][1], spheres[i][2]);
    printf("Distance from origin = %g\n", distance(spheres[i]));
  } 
  delete[] spheres;
}

double ran(){
  return rand()/double(RAND_MAX);
}

double distance(Vector3d v){
  return sqrt((v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
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
