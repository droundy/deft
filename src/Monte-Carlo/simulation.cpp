#include <stdio.h>
#include "Monte-Carlo/monte-carlo.h"

void run(){
  const int N = 13;
  const double R = 1;
  Vector3d *spheres = new Vector3d[N];
  spheres[0] = Vector3d(0,0,0);
  spheres[1] = Vector3d(2*R,0,0);
  spheres[2] = Vector3d(R,R*sqrt(3),0);
  spheres[3] = Vector3d(-R,R*sqrt(3),0);
  spheres[4] = Vector3d(-2*R,0,0);
  spheres[5] = Vector3d(-R,-R*sqrt(3),0);
  spheres[6] = Vector3d(R,-R*sqrt(3),0);
  spheres[7] = Vector3d(R,-R*.5,R*sqrt(3));
  spheres[8] = Vector3d(0,sqrt(3)-.5*R,R*sqrt(3));
  spheres[9] = Vector3d(-R,-.5*R,R*sqrt(3));
  spheres[10] = Vector3d(R,-.5*R,-R*sqrt(3));
  spheres[11] = Vector3d(0,sqrt(3)-.5*R,-R*sqrt(3));
  spheres[12] = Vector3d(-R,-.5*R,-R*sqrt(3)); 
  for(int i = 0; i<N; i++){
    printf("%g  %g  %g\n", spheres[i][0], spheres[i][1], spheres[i][2]);
    printf("Distance from origin = %g\n", distance(spheres[i], Vector3d(0,0,0)));
  }
  FILE *o = fopen("Spheres.dat", "w");
  writeSpheres(spheres, N, o);
   int i = 0;
   int j = 0;
   while(j<100){
   Vector3d temp = move(spheres[i%N], 3, 3, 3);
   if(i%20==0){
      printf("i = %d\n",i);
   }
   if(overlap(spheres, temp, N, R)){
      i++;
  	continue;
  }
  spheres[i] = temp;
  writeSpheres(spheres, N, o);
  j++;
  }   
    spheres[0] = move(spheres[0], 1, 1, 1);
  fclose(o);
  delete[] spheres;
}

bool overlap(Vector3d *spheres, Vector3d v, int n, double R){
  for(int i = 0; i < n; i++){
    if(distance(spheres[n],v)<2*R){
      return true;
    }
  }
    return false;
}

Vector3d move(Vector3d v, double x, double y, double z){
  const double scale =.05;
  Vector3d temp;
  while(true){
    temp = ran3()*scale;
    if((fabs((v+temp)[0]) <= x) && (fabs((v+temp)[1]) <= y) && (fabs((v+temp)[2]) <= z)){
	break;
   }
  }
  // return Vector3d(1,1,1);
  return temp + v;
}
