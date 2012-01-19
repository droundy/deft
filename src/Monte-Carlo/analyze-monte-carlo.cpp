#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//long shell(Vector3d v, long div, double *radius);


int main(int argc, char *argv[]){
  if (argc != 6) {
    printf("usage:  %s radius spheres iterations uncertainty_goal filename \n", argv[0]);
    return 1;
  }
  

  const char *infilename = argv[5];
  const char *outfilenameOrig = argv[5];
  const double rad = atof(argv[1]);
  const int N = atoi(argv[2]);
  const int iterations = atoi(argv[3]);
  const double uncertainty_goal = atof(argv[4]);
  Vector3d *spheres = new Vector3d[N];
  const double R = 1;
  if (uncertainty_goal < 1e-12 || uncertainty_goal > 1.0) {
    printf("Crazy uncertainty goal:  %s\n", argv[1]);
    return 1;
  }
  // We start with randomly-placed spheres, and then gradually wiggle
  // them around until they are all within bounds and are not
  // overlapping.  We do this by creating an "overlap" value which we
  // constrain to never increase.  Note that this may not work at all
  // for high filling fractions, since we could get stuck in a local
  // minimum.
  for(int i=0; i<N; i++){
    spheres[i]=rad*ran3();
  }
  int i = 0;
  clock_t start = clock();
  int num_to_time = 1000;
  int num_timed = 0;

  for(double numOverLaps=countOverLaps(spheres, N, R, rad); numOverLaps>0;){
    if (num_timed++ > num_to_time) {
      clock_t now = clock();
      printf("took %g seconds per initialising iteration\n",
             (now - double(start))/CLOCKS_PER_SEC/num_to_time);
      num_timed = 0;
      start = now;
    }
    Vector3d old =spheres[i%N];
    spheres[i%N]=move(spheres[i%N]);
    double newOverLaps=countOverLaps(spheres, N, R, rad);
    if(newOverLaps>numOverLaps){
      spheres[i%N]=old;
    }
    else{
      numOverLaps=newOverLaps;
    }
    i++;
    if (i%N == 0) {
      printf("numOverLaps=%g\n",numOverLaps);
      //printf("numOverLaps=%g\r",numOverLaps);
      fflush(stdout);
    }
  }
  printf("\nFound initial state!\n");
//////////////////////////////////////////////////////////////////////////
  FILE *o = fopen(infilename, "w");
  fprintf(o,"Radius=%g\n",rad);
  fprintf(o,"Number of Spheres=%d\n", N);
  i = 0;
  long count = 0, workingmoves = 0;
  
  int div = uncertainty_goal*uncertainty_goal*iterations;
  if (div < 10) div = 10;
  printf("Using %d divisions\n", div);
  fflush(stdout);
  printf("%d\t", div);
  printf("%f\t", uncertainty_goal);
  fflush(stdout);
  
  double *radius = new double[div+1];
  for (int l=0;l<div+1;l++) {
    // make each bin have about the same volume
    radius[l] = rad*(pow(double(l)/div, 1.0/3.0) + 0.1*double(i)*uncertainty_goal)/
      (1 + 0.1*div*uncertainty_goal);
  }
  printf("Got to here! ");
  fflush(stdout);
  
  
  //////////////////////////////////////////////////////////////////////////////////////////////
  for (int m=0; m<4;m++){
  
  int *shells = new int[div];
  for (int l=0; l<div; l++) shells[l] = 0;
  
  double *shellsArea = new double [div];
  for (int l=0; l<div; l++) shellsArea[l]=0;
  printf("Got to here now");
  fflush(stdout);
  
  double *density = new double[div];

 // double x;
 // double y;
 // double z;

  double oShellSmall =R+.001*R;
  double oShellMed =R+.01*R;
  double oShellLarge =R+.05*R;
  double oShellGiant =R+.1*R;
  
  double oShellArray[4] = {oShellSmall,oShellMed,oShellLarge,oShellGiant};

  const char *Sizes[4] = {"Small","Med","Large","Giant"};
  
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
  for(int l=0; l<div; l++){
    conShells[l]=0;
    cenConShells[l]=0;
  }
  printf ("%s",outfilename);
  fflush(stdout);
//////////////////////////////////////////////////////////////////////////

  start = clock();
  num_timed = 0;
  for (long j=0; j<iterations; j++){
    if (num_timed++ > num_to_time) {
      clock_t now = clock();
      printf("took %g seconds per iteration\n",
             (now - double(start))/CLOCKS_PER_SEC/num_to_time);
      num_timed = 0;
      start = now;
      // after the first timing, just time things once per percent (as
      // often as we print the % complete messages)
      if (iterations/100 > num_to_time) num_to_time = iterations/100;
    }
  
    count++;
    i++;
    // only write out the sphere positions after they've all had a
    // chance to move, to save on file size.
    if (j%N == 0) {
      //writeSpheres(spheres, N, o);
      for (int j=0;j<N;j++) {
	      shells[shell(spheres[j], div, radius)]++;
		for (long k=0; k<div; k++) {
		  const double rj = distance(spheres[j],Vector3d(0,0,0));
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
		for(long k=0; k<N; k++){
		  for(long count = 0; count<N; count++){
		    if(k!=count && distance(spheres[count],spheres[k]) <= oShell*2){
		  //  for (int cs=0;cs<4;cs++) {
		      //if(j!=count && distance(spheres[count],spheres[j]) <= oShell[cs]*2){
		      conShells[shell(spheres[k],div,radius)]++;
		      cenConShells[shell((spheres[count]+spheres[k])/2,div,radius)]++;
			//  conShells[cs][shell(spheres[j],div,radius)]++;
		      //  cenConShells[cs][shell((spheres[count]+spheres[j])/2,div,radius)]++;
		      }
		    }
		  }
	}
    
//    if (i%N == 0) writeSpheres(spheres, N, o);
    if(j % (iterations/100)==0){
      printf("%g%% complete...\n",j/(iterations*1.0)*100);
      //printf("%g%% complete...\r",j/(times*1.0)*100);
      fflush(stdout);
    }
    Vector3d temp = move(spheres[i%N]);
    if(overlap(spheres, temp, N, R, rad, i%N)){
      continue;
    }
    spheres[i%N] = temp;
    workingmoves++;
  }
  
  printf ("%s",outfilename);
  fflush(stdout);
  //delete[] vecs;
  for(int i=0; i<div; i++){
    printf("Number of spheres in division %d = %d\n", i+1, shells[i]);
  }
  printf ("%s",outfilename);
  fflush(stdout);
  for(int i=0; i<div; i++){
    double rmax = radius[i+1];
    double rmin = radius[i];
    density[i]=shells[i]/(((4/3.*M_PI*rmax*rmax*rmax)-(4/3.*M_PI*rmin*rmin*rmin)))/(iterations/double(N));
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
  //fclose(cData);
  free (outfilename);
  delete[] shells;
  delete[] density;
  delete[] conShells;
  delete[] conDensity;
  delete[] cenConDensity;
  delete[] cenConShells;

  }
  printf("Total number of attempted moves = %ld\n",count);
  printf("Total number of successful moves = %ld\n",workingmoves);
  printf("Acceptance rate = %g\n", workingmoves/double(count));
  delete[] spheres;
  fclose(o);
}


double countOverLaps(Vector3d *spheres, int n, double R, double rad){
  double num = 0;
  for(int j = 0; j<n; j++){
    if(distance(spheres[j],Vector3d(0,0,0))>rad){
      num+=distance(spheres[j],Vector3d(0,0,0))-rad;
    }
    for(int i = j+1; i < n; i++){
      if(distance(spheres[i],spheres[j])<2*R){
        num+=2*R-distance(spheres[i],spheres[j]);
      }
    }
  }
  return num;
}


bool overlap(Vector3d *spheres, Vector3d v, int n, double R, double rad, int s){
  if(distance(v,Vector3d(0,0,0))>rad){
      return true;
  }
  for(int i = 0; i < n; i++){
    if(i==s){
      continue;
    }
    if(distance(spheres[i],v)<2*R){
      return true;
    }
  }
  return false;
}

bool overlap(Vector3d *spheres, Vector3d v, int n, double R, int s, double x, double y, double z){
  if((fabs((v)[0]) > x) || (fabs((v)[1]) > y) || (fabs((v)[2]) > z)){
      return true;
  }
  for(int i = 0; i < n; i++){
    if(i==s){
     continue;
    }
    if(distance(spheres[i],v)<2*R){
      return true;
    }
  }
  return false;
}


Vector3d move(Vector3d v){
  double scale = .5;
  return v+scale*ran3();
}
 
//To be deleted... cvh 
Vector3d move(Vector3d v, double x, double y, double z){
  const double scale =.5;
  Vector3d temp;
  while(true){
    temp = ran3()*scale;
    if((fabs((v+temp)[0]) <= x) && (fabs((v+temp)[1]) <= y) && (fabs((v+temp)[2]) <= z)){
	break;
   }
  }
 return temp + v;
}

bool touch(Vector3d *spheres, Vector3d v, int n, double R, double delta, int s){
  for(int i = 0; i < n; i++){
    if(i==s){
      continue;
    }
    if(distance(spheres[i],v)>=R && distance(spheres[n],v)<=R+delta){
      return true;
    }
  }
  return false;
}

int shell(Vector3d v, int div, double *radius){
  double temp = distance(v,Vector3d(0,0,0));
  for(long count = 0; count<div; count++){
    if(temp<radius[count+1]){
      return count;
    }
  }
  return div-1;
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
