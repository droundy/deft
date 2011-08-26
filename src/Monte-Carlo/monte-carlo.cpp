#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Monte-Carlo/monte-carlo.h"
#include "MersenneTwister.h"

int main(){
  run();
  //  printf("Hello World\n");
  // const int N = 14;
  //Vector3d *spheres = new Vector3d[N];
  //for(int i = 0; i < N; i++){
    //  for(int j = 0; j < 3; j++){
    //  spheres[i][j] = ran();
    // }
  // spheres[i] = ran3();
  //}
  // for(int i = 0; i<N; i++){
  // printf("%g  %g  %g\n", spheres[i][0], spheres[i][1], spheres[i][2]);
  // printf("Distance from origin = %g\n", distance(spheres[i]));
  //} 
  // FILE *o = fopen("Spheres.dat", "w");
  // writeSpheres(spheres, N, o);
  //Need to be able to send in spheres
  //spheres[0] = move(spheres[0], 1, 1, 1);
  //writeSpheres(spheres, N, o);
  //fclose(o);
  //delete[] spheres;
}

double ran(){
  //<<<<<<< HEAD
  MTRand random;
  random.seed();
  return random.rand();
  // return rand()/double(RAND_MAX);
  //=======
  static MTRand my_mtrand;
  return my_mtrand.randExc(); // which is the range of [0,1)
  //>>>>>>> 2511a06ab322dcea9f3f70080c443d0938562932
}

double distance(Vector3d v1, Vector3d v2){
  return sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]));
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
