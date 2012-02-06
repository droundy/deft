#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
using std::vector;

bool periodic_x = true; // will go from -lenx/2 to +lenx/2 
bool periodic_y = true;
bool periodic_z = true;
bool has_x_wall = false;
bool has_y_wall = false;
bool has_z_wall = false;
bool spherical_outer_wall = false;  //spherical walls on outside of entire volume
bool spherical_inner_wall = false;  //sphere at center of entire volume that is a wall
double lenx = 20; 
double leny = 20; 
double lenz = 20; 
double rad = 10;  //of outer spherical walls
double innerRad = 3;  //of inner spherical "solute"
double R = 1;
double oShellSmall =R+.001*R;
double oShellMed =R+.01*R;
double oShellLarge =R+.05*R;
double oShellGiant =R+.1*R;
double oShellArray[4] = {oShellSmall,oShellMed,oShellLarge,oShellGiant};
bool flat_div = false; //the divisions will be equal and will divide from z wall to z wall
Vector3d latx = Vector3d(lenx,0,0);
Vector3d laty = Vector3d(0,leny,0);
Vector3d latz = Vector3d(0,0,lenz);
Vector3d lat[3] = {latx,laty,latz};
bool periodic[3] = {periodic_x, periodic_y, periodic_z};

int main(int argc, char *argv[]){
  if (argc < 5) {
    printf("usage:  %s spheres iterations uncertainty_goal filename \n there will be more!\n", argv[0]);
    return 1;
  }
  
  for (int a=5; a<argc; a+=2){
    printf("Checking a = %d which is %s\n", a, argv[a]);
    if (strcmp(argv[a],"outerSphere") == 0) {
      spherical_outer_wall = true;
      periodic_x = periodic_y = periodic_z = false;
      rad = atof(argv[a+1]);
      printf("Using outerSphere of %g\n", rad);
    }
    if (strcmp(argv[a],"innerSphere") == 0) {
      spherical_inner_wall = true;
      innerRad = atof(argv[a+1]);
    }
    if (strcmp(argv[a],"periodxy") == 0) {
      periodic_x = true;
      periodic_y = true;
      lenx = atof(argv[a+1]);
      leny = atof(argv[a+1]);
    }
    if (strcmp(argv[a],"periodx") == 0) {
      periodic_x = true;
      lenx = atof(argv[a+1]);
    }
    if (strcmp(argv[a],"periody") == 0) {
      periodic_y = true;
      leny = atof(argv[a+1]);
    }
    if (strcmp(argv[a],"periodz") == 0) {
      periodic_z = true;
      lenz = atof(argv[a+1]);
    }
    if (strcmp(argv[a],"wallx") == 0) {
      has_x_wall = true;
      lenx = atof(argv[a+1]);
      periodic_x = false;
    }
    if (strcmp(argv[a],"wally") == 0) {
      has_y_wall = true;
      lenx = atof(argv[a+1]);
      periodic_y = false;
    }
    if (strcmp(argv[a],"wallz") == 0) {
      has_z_wall = true;
      lenx = atof(argv[a+1]);
      periodic_z = false;
    }
    if (strcmp(argv[a],"flatdiv") == 0) {
      flat_div = true; //otherwise will default to radial divisions
      a -= 1;
    }
  }
  printf("flatdiv = %s\n", flat_div ? "true" : "false");
  
  const char *outfilename = argv[4];
  printf ("this is %s",outfilename);
  fflush(stdout);
  const int N = atoi(argv[1]);
  const int iterations = atoi(argv[2]);
  const double uncertainty_goal = atof(argv[3]);
  Vector3d *spheres = new Vector3d[N];
  if (uncertainty_goal < 1e-12 || uncertainty_goal > 1.0) {
    printf("Crazy uncertainty goal:  %s\n", argv[1]);
    return 1;
  }

  //////////////////////////////////////////////////////////////////////////////////////////
  // We start with randomly-placed spheres, and then gradually wiggle
  // them around until they are all within bounds and are not
  // overlapping.  We do this by creating an "overlap" value which we
  // constrain to never increase.  Note that this may not work at all
  // for high filling fractions, since we could get stuck in a local
  // minimum.
  for(int i=0; i<N; i++){
    spheres[i]=rad*ran3();
  }
  clock_t start = clock();
  int num_to_time = 1000;
  int num_timed = 0;
  int i = 0;
  double scale = .5;
  
  for(double numOverLaps=countOverLaps(spheres, N, R); numOverLaps>0;){
    if (num_timed++ > num_to_time) {
      clock_t now = clock();
      printf("took %g seconds per initialising iteration\n",
             (now - double(start))/CLOCKS_PER_SEC/num_to_time);
      num_timed = 0;
      start = now;
    }
    Vector3d old =spheres[i%N];
    spheres[i%N]=move(spheres[i%N],scale);
    double newOverLaps=countOverLaps(spheres, N, R);
    if(newOverLaps>numOverLaps){
      spheres[i%N]=old;
      scale = scale*0.98;
    } else if (scale < 5) {
      numOverLaps=newOverLaps;
      scale = scale*1.02;
    }
    i++;
    if (i%N == 0) {
      printf("numOverLaps=%g\n",numOverLaps);
      fflush(stdout);
    }
  }
  printf("\nFound initial state!\n");
//////////////////////////////////////////////////////////////////////////
  //FILE *o = fopen(outfilename, "w");
  //fprintf(o,"Radius=%g\n",rad);
  //fprintf(o,"Number of Spheres=%d\n", N);
  
  int div = uncertainty_goal*uncertainty_goal*iterations;
  if (div < 10) div = 10;
  printf("Using %d divisions\n", div);
  fflush(stdout);
  printf("%d\t", div);
  printf("%f\t", uncertainty_goal);
  fflush(stdout);
  
  double *radius = new double[div+1];
  double *sections = new double [div+1];

  if (flat_div){
    double size = lenz/div;
    for (int s=0; s<div+1; s++){
      sections[s] = size*s - lenz/2.0;
    }
  } else {
    if (spherical_inner_wall){
      vector <double> Radius;
      int newdiv = div;
      int rs = 0;
      for (int l=0;l<div+1;l++) {
	Radius.push_back(rad*(pow(double(l)/(div), 1.0/3.0) + 0.1*double(l)*uncertainty_goal)/
			 (1 + 0.1*(div)*uncertainty_goal) );
      }
      for (int ks=0;ks<1000000;ks++){
	while (Radius[rs] < innerRad) rs++;
	newdiv = newdiv+(rs-1);
	Radius.clear();
	for (int l=0;l<newdiv+1;l++) {
	  Radius.push_back (rad*(pow(double(l)/(newdiv), 1.0/3.0) + 0.1*double(l)*uncertainty_goal)/
			    (1 + 0.1*(newdiv)*uncertainty_goal) );
	}
	if (newdiv-(rs-1) == div) break;
	printf("HEEEEEEEEEEERE    %d",ks);
	fflush(stdout);
      }
      rs=0;
      while (Radius[rs]<innerRad)rs++;
      double VolDelta = pow(Radius[rs],3)-pow(innerRad,3);
      double VolDeltaNext = pow(Radius[rs+1],3)-pow(Radius[rs],3);
      if (VolDelta < VolDeltaNext/2.0){
	Radius[rs] = pow( (pow(Radius[rs+1],3)-pow(innerRad,3))/2 + pow(innerRad,3) ,1.0/3.0);
      }
      radius[0] = innerRad;
      for (int k=0;k<div;k++){
	radius[k+1]= Radius[k+rs];
      }
    } else {
      for (int l=0;l<div+1;l++) {
	// make each bin have about the same volume
	radius[l] = rad*(pow(double(l)/(div), 1.0/3.0) + 0.1*double(l)*uncertainty_goal)/
	  (1 + 0.1*(div)*uncertainty_goal);
      }
    }
  }
  
  if (flat_div){
    for (int k=0;k<div+1;k++) printf("section  = %f\n",sections[k]);
    fflush(stdout);
  } else{
    for (int k=0;k<div+1;k++) printf("radius  = %f\n",radius[k]);
    fflush(stdout);
  }


  printf("Look HEREE WAAAAAAAAAAAAAA !!!!%f\n",radius[0]); 
  printf("%f\n",radius[1]); 
  printf("%f\n",radius[2]);
  fflush(stdout); 
  //////////////////////////////////////////////////////////////////////////////////////////////
  scale = .5;
  int count = 0;
  int *shells = new int[div];
  for (int l=0; l<div; l++) shells[l] = 0;
  
  double *shellsArea = new double [div];
  for (int l=0; l<div; l++) shellsArea[l]=0;
  
  double *density = new double[div];
  
  double *SconDensity = new double[div]; double *ScenConDensity = new double[div];
  int *SconShells = new int[div]; int *ScenConShells = new int[div];
  double *MconDensity = new double[div]; double *McenConDensity = new double[div];
  int *MconShells = new int[div]; int *McenConShells = new int[div];
  double *LconDensity = new double[div]; double *LcenConDensity = new double[div];
  int *LconShells = new int[div]; int *LcenConShells = new int[div];
  double *GconDensity = new double[div]; double *GcenConDensity = new double[div];
  int *GconShells = new int[div]; int *GcenConShells = new int[div];
  
  
  for(int l=0; l<div; l++){
    SconShells[l]=0; MconShells[l]=0; LconShells[l]=0; GconShells[l]=0;
    ScenConShells[l]=0;McenConShells[l]=0;LcenConShells[l]=0;GcenConShells[l]=0;
  }
  /////////////////////////////////////////////////////////////////////////////
  
  start = clock();
  num_timed = 0;
  int workingmoves=0;
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
    // only write out the sphere positions after they've all had a
    // chance to move
    if (j%N == 0) {
      for (int i=0;i<N;i++) {
	shells[shell(spheres[i], div, radius, sections)]++;
	if (!flat_div){
	  for (long k=0; k<div; k++) {
	    const double ri = distance(spheres[i],Vector3d(0,0,0));
	    if (ri < radius[k+1] + R && ri + radius[k+1] > R && ri > radius[k] - R) {
	      // There is at least some overlap with shell k! (not so easy)
	      double costhetamax, costhetamin;
	      if (ri > radius[k] + R) {
		costhetamin = 1;
	      } else if (radius[k] + ri < R) {
		costhetamin = 1;
	      } else {
		costhetamin = (ri*ri - radius[k]*radius[k] + R*R)/(2*ri*R);
	      }
	      if (ri < radius[k+1] - R) {
		costhetamax = -1;
	      } else {
		costhetamax = (ri*ri - radius[k+1]*radius[k+1] + R*R)/(2*ri*R);
	      }
	      assert(costhetamin >= costhetamax);
	      shellsArea[k] += 2*M_PI*R*R*(costhetamin-costhetamax);
	    }
	  }
	} else {
	  for (int k=0; k<div+1; k++){
	    double dl = spheres[i][2] - sections[k];
	    double dh = spheres[i][2] - sections[k+1];
	    if (dl > R) dl = R;
	    if (dl < -R) dl = -R;
	    if (dh > R) dh = R;
	    if (dh < -R) dh = -R;
	    shellsArea[k] += 2*M_PI*R*(dl-dh);
	  }
	}
      }	
      for(long k=0; k<N; k++){
	for(long n = 0; n<N; n++){
	  //if(k!=n && distance(spheres[n],spheres[k]) <= oShell*2 || cornTouch(spheres[n],spheres[k],R)){
	  if(k!=n && touch(spheres[n],spheres[k],oShellArray[0])){
	    SconShells[shell(spheres[k],div,radius,sections)]++;
	    ScenConShells[shell((spheres[n]+spheres[k])/2,div,radius,sections)]++;
	  } 
	  if(k!=n && touch(spheres[n],spheres[k],oShellArray[1])){
	    MconShells[shell(spheres[k],div,radius,sections)]++;
	    McenConShells[shell((spheres[n]+spheres[k])/2,div,radius,sections)]++;
	  }
	  if(k!=n && touch(spheres[n],spheres[k],oShellArray[2])){
	    LconShells[shell(spheres[k],div,radius,sections)]++;
	    LcenConShells[shell((spheres[n]+spheres[k])/2,div,radius,sections)]++;
	  }
	  if(k!=n && touch(spheres[n],spheres[k],oShellArray[3])){
	    GconShells[shell(spheres[k],div,radius,sections)]++;
	    GcenConShells[shell((spheres[n]+spheres[k])/2,div,radius,sections)]++;
	  }
	}
      }
    }
    if(j % (iterations/100)==0){
      printf("%g%% complete...\n",j/(iterations*1.0)*100);
      //printf("%g%% complete...\r",j/(times*1.0)*100);
      fflush(stdout);
    }
    Vector3d temp = move(spheres[j%N],scale);
    count++;
    if(overlap(spheres, temp, N, R, j%N)){
      scale = scale*(0.98);
      continue;
    }
    spheres[j%N] = temp;
    workingmoves++;
    if (scale < 5) scale = scale*1.02;
  }

//////////////////////////////////////////////////////////////////////////////////////////
  
  for(int i=0; i<div; i++){
    printf("Number of spheres in division %d = %d\n", i+1, shells[i]);
  }
  if (!flat_div){
    for(int i=0; i<div; i++){
      double rmax = radius[i+1];
      double rmin = radius[i];
      density[i]=shells[i]/(((4/3.*M_PI*rmax*rmax*rmax)-(4/3.*M_PI*rmin*rmin*rmin)))/(iterations/double(N));
    }
  } else {
    for(int i=0; i<div; i++){
      density[i]=shells[i]/(lenx*leny*lenz/div)/(iterations/double(N));
    }
  }
  
  for(int i=0; i<div; i++){
    SconDensity[i]=((SconShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShellArray[0]*8*oShellArray[0]*oShellArray[0]-4/3.*M_PI*8*R*R*R));
    ScenConDensity[i]=4*M_PI*R*R*((ScenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShellArray[0]*oShellArray[0]*oShellArray[0]-4/3.*M_PI*8*R*R*R));
    MconDensity[i]=((MconShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShellArray[1]*8*oShellArray[1]*oShellArray[1]-4/3.*M_PI*8*R*R*R));
    McenConDensity[i]=4*M_PI*R*R*((McenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShellArray[1]*oShellArray[1]*oShellArray[1]-4/3.*M_PI*8*R*R*R));
    LconDensity[i]=((LconShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShellArray[2]*8*oShellArray[2]*oShellArray[2]-4/3.*M_PI*8*R*R*R));
    LcenConDensity[i]=4*M_PI*R*R*((LcenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShellArray[2]*oShellArray[2]*oShellArray[2]-4/3.*M_PI*8*R*R*R));
    GconDensity[i]=((GconShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShellArray[3]*8*oShellArray[3]*oShellArray[3]-4/3.*M_PI*8*R*R*R));
    GcenConDensity[i]=4*M_PI*R*R*((GcenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShellArray[3]*oShellArray[3]*oShellArray[3]-4/3.*M_PI*8*R*R*R));
  }
  for(int i=0; i<div; i++){
    printf("Number of contacts in division %d = %d%d%d%d\n", i+1, 
	   SconShells[i], MconShells[i], LconShells[i], GconShells[i]);
    printf("Number of contacts (center) in division %d = %d%d%d%d\n", i+1, 
	   ScenConShells[i], McenConShells[i], LcenConShells[i], GcenConShells[i]);
  }
  
//FILE *out = fopen((const char *)outfilename,"w");
  FILE *out = fopen((const char *)outfilename,"w");
  if (out == NULL) {
    printf("Error creating file %s\n", outfilename);
    return 1;
  }
  if (flat_div){
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 0.5*(sections[0]+sections[1]), density[0],
	    SconDensity[0], ScenConDensity[0], MconDensity[0], McenConDensity[0], LconDensity[0], LcenConDensity[0], GconDensity[0], GcenConDensity[0]);
  } else if (spherical_inner_wall) {
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", radius[0], 0.0, 
	    SconDensity[0], ScenConDensity[0], MconDensity[0], McenConDensity[0], LconDensity[0], LcenConDensity[0], GconDensity[0], GcenConDensity[0]);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 0.5*(radius[0]+radius[1]), density[0],
	    SconDensity[0], ScenConDensity[0], MconDensity[0], McenConDensity[0], LconDensity[0], LcenConDensity[0], GconDensity[0], GcenConDensity[0]);
  } else {
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" , 0.0, density[0], 
	    SconDensity[0], ScenConDensity[0], MconDensity[0], McenConDensity[0], LconDensity[0], LcenConDensity[0], GconDensity[0], GcenConDensity[0]);
  }
  
  int divtoprint = div;
  if (!spherical_outer_wall) divtoprint = div - 1;
  if (!flat_div) {
    for(int i=1; i<divtoprint; i++){
      fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	      0.5*(radius[i]+radius[i+1]), density[i], 
	      SconDensity[i], ScenConDensity[i], MconDensity[i], McenConDensity[i], LconDensity[i], LcenConDensity[i], GconDensity[i], GcenConDensity[i]);
    }
  } else {
    for(int i=1; i<div; i++){
      fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
	      0.5*(sections[i]+sections[i+1]), density[i], 
	      SconDensity[i], ScenConDensity[i], MconDensity[i], McenConDensity[i], LconDensity[i], LcenConDensity[i], GconDensity[i], GcenConDensity[i]);
    }
  }

  printf("Total number of attempted moves = %d\n",count);
  printf("Total number of successful moves = %d\n",workingmoves);
  printf("Acceptance rate = %g\n", workingmoves/double(count));
  //delete[] shells;
  printf("herefirst");
  fflush(stdout);
  //delete[] density;
  printf("heresecond");
  fflush(stdout);
  delete[] SconShells; 
  printf("herethird");
  fflush(stdout);
  delete[] SconDensity; 
  printf("herefourth");
  fflush(stdout);
  delete[] ScenConDensity; delete[] ScenConShells;
  delete[] MconShells; delete[] MconDensity; delete[] McenConDensity; delete[] McenConShells;
  delete[] LconShells; delete[] LconDensity; delete[] LcenConDensity; delete[] LcenConShells;
  delete[] GconShells; delete[] GconDensity; delete[] GcenConDensity; delete[] GcenConShells;
  printf("herefifth");
  fflush(stdout);
  delete[] spheres;
  printf("heresixth\n");
  fflush(stdout);
  fclose(out);
}


double countOverLaps(Vector3d *spheres, int n, double R){
  double num = 0;
  for(int j = 0; j<n; j++){
    for(int i = j+1; i < n; i++){
        if(distance(spheres[i],spheres[j])<2*R){
          num+=2*R-distance(spheres[i],spheres[j]);
	}
    }
    if (spherical_outer_wall){
      if(distance(spheres[j],Vector3d(0,0,0))>rad){
        num += distance(spheres[j],Vector3d(0,0,0))-rad;
      }
    }
    if (spherical_inner_wall){
      if(distance(spheres[j],Vector3d(0,0,0))<innerRad){
        num -= distance(spheres[j],Vector3d(0,0,0))-innerRad;
      }
    }
    if (has_x_wall && periodic_x){
      if (spheres[j][0] > lenx/2 ){
	num += spheres[j][0]-(lenx/2);
      } else if (spheres[j][0] < -lenx/2){
	num -= spheres[j][0] + (lenx/2);
      }
    }
    if (has_y_wall && periodic_y){
      if (spheres[j][1] > leny/2 ){
	num += spheres[j][1]-(leny/2);
      } else if (spheres[j][1] < -leny/2){
	num -= spheres[j][1] + (leny/2);
      }
    }
    if (has_z_wall && periodic_z){
      if (spheres[j][2] > lenz/2 ){
	num += spheres[j][2]-(lenz/2);
      } else if (spheres[j][2] < -lenz/2){
	num -= spheres[j][2] + (lenz/2);
      }
    }
    Vector3d lat[3] = {latx,laty,latz};
    bool periodic[3] = {periodic_x, periodic_y, periodic_z};
    for(int i = j+1; i < n; i++){
      for (int k=0; k<3; k++){
	if (periodic[k]){
	  if (distance(spheres[j],spheres[i]+lat[k]) < 2*R){
	    num += 2*R - distance(spheres[j],spheres[i]+lat[k]);
	  } else if (distance(spheres[j],spheres[i]-lat[k]) < 2*R) {
	    num += 2*R - distance(spheres[j],spheres[i]-lat[k]);
	  }
	} 
	for (int m=k+1; m<3; m++){
	  if (periodic[m] && periodic[k]){
	    if (distance(spheres[j],spheres[i]+lat[k]+lat[m]) < 2*R){
	      num += 2*R - distance(spheres[j],spheres[i]+latx+laty);
	    }
	    if (distance(spheres[j],spheres[i]-lat[k]-lat[m]) < 2*R){
	      num += 2*R - distance(spheres[j],spheres[i]+latx+laty);
	    }
	    if (distance(spheres[j],spheres[i]+lat[k]-lat[m]) < 2*R){
	      num += 2*R - distance(spheres[j],spheres[i]+latx+laty);
	    }
	    if (distance(spheres[j],spheres[i]-lat[k]+lat[m]) < 2*R){
	      num += 2*R - distance(spheres[j],spheres[i]+latx+laty);
	    }
	  }
	}
      }
      if (periodic[0] && periodic[1] && periodic[2]){
	if (distance(spheres[j],spheres[i]+latx+laty+latz) < 2*R){
	  num += 2*R - distance(spheres[j],spheres[i]+latx+laty+latz);
	}
	if (distance(spheres[j],spheres[i]+latx+laty-latz) < 2*R){
	  num += 2*R - distance(spheres[j],spheres[i]+latx+laty-latz);
	}
	if (distance(spheres[j],spheres[i]+latx-laty+latz) < 2*R){
	  num += 2*R - distance(spheres[j],spheres[i]+latx-laty+latz);
	}
	if (distance(spheres[j],spheres[i]-latx+laty+latz) < 2*R){
	  num += 2*R - distance(spheres[j],spheres[i]-latx+laty+latz);
	}
	if (distance(spheres[j],spheres[i]-latx-laty+latz) < 2*R){
	  num += 2*R - distance(spheres[j],spheres[i]-latx-laty+latz);
	}
	if (distance(spheres[j],spheres[i]-latx+laty-latz) < 2*R){
	  num += 2*R - distance(spheres[j],spheres[i]-latx+laty-latz);
	}
	if (distance(spheres[j],spheres[i]+latx-laty-latz) < 2*R){
	  num += 2*R - distance(spheres[j],spheres[i]+latx-laty-latz);
	}
	if (distance(spheres[j],spheres[i]-latx-laty-latz) < 2*R){
	  num += 2*R - distance(spheres[j],spheres[i]-latx-laty-latz);
	}
      }
    }
  }
  return num;
}


bool overlap(Vector3d *spheres, Vector3d v, int n, double R, int s){
  if (spherical_outer_wall){
    if (distance(v,Vector3d(0,0,0)) > rad) return true;
  }
  if (spherical_inner_wall) {
    if (distance(v,Vector3d(0,0,0)) < innerRad) return true;
  }
  if (has_x_wall){
    if (v[0] > lenx/2 || v[0] < -lenx/2) return true;
  }
  if (has_y_wall){
    if (v[1] > leny/2 || v[1] < -leny/2) return true;
  }
  if (has_z_wall){
    if (v[2] > lenz/2 || v[2] < -lenz/2) return true;
  }
  for(int i = 0; i < n; i++){
    if(i!=s){
      if(distance(spheres[i],v)<2*R){
	return true;
      }
      for (int k=0; k<3; k++){
	if (periodic[k]){
	  if (distance(v,spheres[i]+lat[k]) < 2*R) return true;
	  if (distance(v,spheres[i]-lat[k]) < 2*R) return true;
	} 
	for (int m=k+1; m<3; m++){
	  if (periodic[m] && periodic[k]){
	    if (distance(v,spheres[i]+lat[k]+lat[m]) < 2*R) return true;
	    if (distance(v,spheres[i]-lat[k]-lat[m]) < 2*R) return true;
	    if (distance(v,spheres[i]+lat[k]-lat[m]) < 2*R) return true;
	    if (distance(v,spheres[i]-lat[k]+lat[m]) < 2*R) return true;
	  }
	}
      }
      if (periodic[0] && periodic[1] && periodic[2]){
	if (distance(v,spheres[i]+latx+laty+latz) < 2*R) return true;
	if (distance(v,spheres[i]+latx+laty-latz) < 2*R) return true;
	if (distance(v,spheres[i]+latx-laty+latz) < 2*R) return true;
	if (distance(v,spheres[i]-latx+laty+latz) < 2*R) return true;
	if (distance(v,spheres[i]-latx-laty+latz) < 2*R) return true;
	if (distance(v,spheres[i]-latx+laty-latz) < 2*R) return true;
	if (distance(v,spheres[i]+latx-laty-latz) < 2*R) return true;
	if (distance(v,spheres[i]-latx-laty-latz) < 2*R) return true;
      }
    }
  }
  return false;
}


Vector3d move(Vector3d v,double scale){
  Vector3d newv = v+scale*ran3();  
  if (periodic_x){
    while (newv[0] > lenx/2){
      newv[0] -= lenx;
    }
    while (newv[0] < -lenx/2){
      newv[0] += lenx;
    }
  }
  if (periodic_y){
    while (newv[1] > leny/2){
      newv[1] -= leny;
    }
    while (newv[1] < -leny/2){
      newv[1] += leny;
    }
  }
  if (periodic_z){
    while (newv[2] > lenz/2){
      newv[2] -= lenz;
    }
    while (newv[2] < -lenz/2){
      newv[2] += lenz;
    }
  }
  return newv;
}
 

bool touch(Vector3d w, Vector3d v, double oShell){
  if (distance(v,w) < 2*oShell) return true;
  for (int k=0; k<3; k++){
	if (periodic[k]){
	  if (distance(v,w+lat[k]) < 2*oShell) return true;
	  if (distance(v,w-lat[k]) < 2*oShell) return true;
	} 
	for (int m=k+1; m<3; m++){
	  if (periodic[m] && periodic[k]){
	    if (distance(v,w+lat[k]+lat[m]) < 2*oShell) return true;
	    if (distance(v,w-lat[k]-lat[m]) < 2*oShell) return true;
	    if (distance(v,w+lat[k]-lat[m]) < 2*oShell) return true;
	    if (distance(v,w-lat[k]+lat[m]) < 2*oShell) return true;
	  }
	}
  }
      if (periodic[0] && periodic[1] && periodic[2]){
	if (distance(v,w+latx+laty+latz) < 2*oShell) return true;
	if (distance(v,w+latx+laty-latz) < 2*oShell) return true;
	if (distance(v,w+latx-laty+latz) < 2*oShell) return true;
	if (distance(v,w-latx+laty+latz) < 2*oShell) return true;
	if (distance(v,w-latx-laty+latz) < 2*oShell) return true;
	if (distance(v,w-latx+laty-latz) < 2*oShell) return true;
	if (distance(v,w+latx-laty-latz) < 2*oShell) return true;
	if (distance(v,w-latx-laty-latz) < 2*oShell) return true;
      }
      return false;
}



int shell(Vector3d v, int div, double *radius, double *sections){
  if (!flat_div){ 
    double temp = distance(v,Vector3d(0,0,0));
    for(long count = 0; count<div; count++){
      if(temp<radius[count+1]) return count;
    }
  } else {
    double temp = v[2];
    for(long count = 0; count<div; count++){
      if(temp < sections[count+1]) return count;
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
