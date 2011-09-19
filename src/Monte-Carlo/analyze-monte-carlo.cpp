#include <stdio.h>
#include "Monte-Carlo/monte-carlo.h"

int main(){
  const int div = 10;
  const double R = 1;
  int *shells = new int[div];
  double *density = new double[div];

  int n = 0;
  while(n<div){
    shells[n]=0;
    n++;
  }

  double x;
  double y;
  double z;

  FILE *data;
  data = fopen("Spheres(120-in-6).dat", "r");
  double rad;
  int N;
  if(fscanf(data, " Radius = %lg",&rad)!=1){
    printf("Bad file format");
    exit(1);
  }
  if(fscanf(data, " Number of Spheres = %d",&N)!=1){
    printf("Bad file format");
    exit(1);
  }
  int iterations = 0;
  while(fscanf(data, "%lg %lg %lg", &x, &y, &z) == 3){
    Vector3d * v = new Vector3d(x,y,z);
    shells[shell(*v,div,rad)]++;
    iterations++;
  }
  fclose(data);

  double *conDensity = new double[div];
  double *cenConDensity = new double[div];
  int *conShells = new int[div];
  int *cenConShells = new int[div];
  for(int i=0; i<div; i++){
    conShells[i]=0;
    cenConShells[i]=0;
  }
  double oShell =R+.1*R;
  FILE *cData;
  cData = fopen("Spheres(120-in-6).dat","r");

  if(fscanf(cData, " Radius = %lg",&rad)!=1){
    printf("Bad file format");
    exit(1);
  }
  if(fscanf(cData, " Number of Spheres = %d",&N)!=1){
    printf("Bad file format");
    exit(1);
  }

  Vector3d *vecs = new Vector3d[N];
  for(int i=0; i<(iterations/N); i++){
    for(int j=0; j<N; j++){
      if(fscanf(cData, "%lg %lg %lg", &x, &y, &z)!=3){
	printf("Bad data vector\n");
	exit(1);
      }
      vecs[j] = Vector3d(x,y,z);
    }
    for(int j=0; j<N; j++){
      for(int count = 0; count<N; count++){
	if(j!=count && distance(vecs[count],vecs[j]) <= oShell*2){
	  conShells[shell(vecs[j],div,rad)]++;
	  cenConShells[shell((vecs[count]+vecs[j])/2,div,rad)]++;
	}
      }
    }
  }
  for(int i=0; i<div; i++){
    printf("Number of spheres in division %d = %d\n", i+1, shells[i]);
  }
  for(int i=0; i<div; i++){
    double rmax = ((i+1)*rad/div);
    double rmin = (i*rad/div);
    density[i]=shells[i]/(((4/3.*M_PI*rmax*rmax*rmax)-(4/3.*M_PI*rmin*rmin*rmin)))/(iterations/N);
  }
  for(int i=0; i<div; i++){
    conDensity[i]=((conShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShell*8*oShell*oShell-4/3.*M_PI*8*R*R*R));
    cenConDensity[i]=((cenConShells[i]+0.0)/shells[i])/((4/3.*M_PI*8*oShell*oShell*oShell-4/3.*M_PI*8*R*R*R));
  }
  for(int i=0; i<div; i++){
    printf("Number of contacts in division %d = %d\n", i+1, conShells[i]);
    printf("Number of contacts (center) in division %d = %d\n", i+1, cenConShells[i]);
  }
  for(int i=0; i<div; i++){
    printf("Contact density in division %d = %g \n", i+1, conDensity[i]);
    printf("Contact density (center) in division %d = %g \n", i+1, cenConDensity[i]);
    printf("Density in division %d = %g \n", i+1, density[i]);
    printf("Filling fraction in division %d = %g \n", i+1, density[i]*(4/3.*M_PI));
  }
  FILE *out = fopen("density.dat","w");
  for(int i=0; i<div; i++){
    fprintf(out, "%g\t%g\t%g\t%g\n",(i+.5)*rad/div, density[i], conDensity[i], cenConDensity[i]);
  }
  fclose(out);
  fclose(cData);
  delete[] shells;
  delete[] density;
  delete[] conShells;
  delete[] conDensity;
  delete[] cenConDensity;
  delete[] cenConShells;
  return 0;
}

int shell(Vector3d v, int div, double rad){
  double temp = distance(v,Vector3d(0,0,0));
  int count = 1;
  for(double i = rad/div; i <= div; i+=rad/div){
    if(count==1){
      if(temp <= i){
	return count;
      }
      count++;
      continue;
    }
    if(temp > i && temp <= (i+rad/div)){
      return count;
    }
    count++;
  }
  return 0;
}

double distance(Vector3d v1, Vector3d v2){
  return sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]));
}

