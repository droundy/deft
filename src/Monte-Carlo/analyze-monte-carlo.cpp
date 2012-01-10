//Changed?
#include <stdio.h>
#include "Monte-Carlo/monte-carlo.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]){
  if (argc != 4) {
    printf("usage:  %s uncertainty_goal inputfilename outputfilename\n", argv[0]);
    return 1;
  }
  const char *infilename = argv[2];
  const char *outfilenameOrig = argv[3];
  const double uncertainty_goal = atof(argv[1]);
  printf("Analyzing file %s with %g precision...\n", infilename, uncertainty_goal);
  
  if (uncertainty_goal < 1e-12 || uncertainty_goal > 1.0) {
    printf("Crazy uncertainty goal:  %s\n", argv[1]);
    return 1;
  }
  const double R = 1;
   
  // First, we'll open the file and count the number of iterations
  FILE *data;
  data = fopen(infilename, "r");
  if (data == NULL) {
    printf("Error opening file %s!\n", infilename);
    return 1;
  }
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
  {
    double dummy;
    while(fscanf(data, "%lg %lg %lg", &dummy, &dummy, &dummy) == 3) iterations++;
  }
  fclose(data);

  int div = pow(uncertainty_goal*uncertainty_goal*iterations, 1.0/3.0);
  if (div < 10) div = 10;
  printf("Using %d divisions\n", div);
  fflush(stdout);

  double *radius = new double[div+1];
  for (int i=0;i<div+1;i++) {
    // make each bin have about the same volume
    radius[i] = rad*(pow(double(i)/div, 1.0/3.0) + 0.1*double(i)*uncertainty_goal)/
      (1 + 0.1*div*uncertainty_goal);
  }
  printf("Got to here! ");
  fflush(stdout);
  
  for (int m=0; m<4;m++){
  
  int *shells = new int[div];
  for (int i=0; i<div; i++) shells[i] = 0;
  
  double *shellsArea = new double [div];
  for (int i=0; i<div; i++) shellsArea[i]=0;
  printf("Got to here now");
  fflush(stdout);
  double *density = new double[div];

  double x;
  double y;
  double z;

  double oShellSmall =R+.1*R;
  double oShellMed =R+.1*R;
  double oShellLarge =R+.1*R;
  double oShellGiant =R+.1*R;
  
  double oShellArray[4] = {oShellSmall,oShellMed,oShellLarge,oShellGiant};

  const char *Sizes[4] = {"Small","Med","Large","Giant"};

  FILE *cData;
  cData = fopen(infilename,"r");
  
  if(fscanf(cData, " Radius = %lg",&rad)!=1){
    printf("Bad file format");
    exit(1);
  }
  if(fscanf(cData, " Number of Spheres = %d",&N)!=1){
    printf("Bad file format");
    exit(1);
  }
  fflush(stdout);
  
  double oShell = oShellArray[m];
  char *outfilename = (char *) malloc (sizeof(char) * (strlen(outfilenameOrig) + strlen(Sizes[m])));
  strcpy (outfilename,Sizes[m]);
  strcat (outfilename,outfilenameOrig);
  printf ("%s",outfilename);
  fflush(stdout);

  double *conDensity = new double[div];
  double *cenConDensity = new double[div];
  int *conShells = new int[div];
  int *cenConShells = new int[div];
  for(int i=0; i<div; i++){
    conShells[i]=0;
    cenConShells[i]=0;
  }
  printf ("%s",outfilename);
  fflush(stdout);
  
  Vector3d *vecs = new Vector3d[N];
  for(int i=0; i<(iterations/N); i++){
    for(int j=0; j<N; j++){
      if(fscanf(cData, "%lg %lg %lg", &x, &y, &z)!=3){
        printf("Bad data vector\n");
        exit(1);
      }
      vecs[j] = Vector3d(x,y,z);
      shells[shell(vecs[j], div, radius)]++;
      for (int k=0; k<div; k++) {
	const double rj = distance(vecs[j],Vector3d(0,0,0));
	if (rj < radius[k+1] + R && rj + radius[k+1] > R && rj > radius[k] - R) {
	  // There is at least some overlap with shell k! (not so easy)
	  double costhetamax, costhetamin;
	  if (rj > radius[k] + R) {
	    costhetamin = 1;
	  } else if (radius[k] + rj < R) {
	    costhetamin = 1;
	  } else {
	    costhetamin = (rj*rj - radius[k]*radius[k] + R*R)/(2*rj*R);
	  }
	  if (rj < radius[k+1] - R) {
	    costhetamax = -1;
	  } else {
	    costhetamax = (rj*rj - radius[k+1]*radius[k+1] + R*R)/(2*rj*R);
	  }
	  assert(costhetamin >= costhetamax);
	  shellsArea[k] += 2*M_PI*R*R*(costhetamin-costhetamax);
	}
      }
    }
    printf ("%s",outfilename);
    fflush(stdout);
    for(int j=0; j<N; j++){
      for(int count = 0; count<N; count++){
        if(j!=count && distance(vecs[count],vecs[j]) <= oShell*2){
          conShells[shell(vecs[j],div,radius)]++;
          cenConShells[shell((vecs[count]+vecs[j])/2,div,radius)]++;
        }
      }
    }
  }
  printf ("%s",outfilename);
  fflush(stdout);
  delete[] vecs;
  for(int i=0; i<div; i++){
    printf("Number of spheres in division %d = %d\n", i+1, shells[i]);
  }
  printf ("%s",outfilename);
  fflush(stdout);
  for(int i=0; i<div; i++){
    double rmax = radius[i+1];
    double rmin = radius[i];
    density[i]=shells[i]/(((4/3.*M_PI*rmax*rmax*rmax)-(4/3.*M_PI*rmin*rmin*rmin)))/(iterations/N);
  }
  printf ("%s",outfilename);
  fflush(stdout);
  for(int i=0; i<div; i++){
    conDensity[i]=((conShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShell*8*oShell*oShell-4/3.*M_PI*8*R*R*R));
    cenConDensity[i]=4*M_PI*R*R*((cenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShell*oShell*oShell-4/3.*M_PI*8*R*R*R));
  }
  for(int i=0; i<div; i++){
    printf ("%s",outfilename);
    fflush(stdout);
    printf("Number of contacts in division %d = %d\n", i+1, conShells[i]);
    printf("Number of contacts (center) in division %d = %d\n", i+1, cenConShells[i]);
  }
  printf ("%s",outfilename);
  FILE *out = fopen((const char *)outfilename,"w");
  printf ("%s",outfilename);
  if (out == NULL) {
    printf("Error creating file %s\n", outfilename);
    return 1;
  }
  fprintf(out, "%g\t%g\t%g\t%g\n", 0.0, density[0], conDensity[0], cenConDensity[0]);
  for(int i=1; i<div; i++){
    fprintf(out, "%g\t%g\t%g\t%g\n",
            0.5*(radius[i]+radius[i+1]), density[i], conDensity[i], cenConDensity[i]);
  }
  fclose(out);
  fclose(cData);
  free (outfilename);
  delete[] shells;
  delete[] density;
  delete[] conShells;
  delete[] conDensity;
  delete[] cenConDensity;
  delete[] cenConShells;
  }
  return 0;
}

int shell(Vector3d v, int div, double *radius){
  double temp = distance(v,Vector3d(0,0,0));
  for(int count = 0; count<div; count++){
    if(temp<radius[count+1]){
      return count;
    }
  }
  return div-1;
}
