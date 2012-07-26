#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
using std::vector;

long shell(Vector3d v, long div, double *radius, double *sections);
double countOverLaps(Vector3d *spheres, long n, double R);
double countOneOverLap(Vector3d *spheres, long n, long j, double R);
bool overlap(Vector3d *spheres, Vector3d v, long n, double R, long s);

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
double oShellSmall =R+.0001*R;
double oShellMed =R+.0005*R;
double oShellLarge =R+.001*R;
double oShellGiant =R+.01*R;
double oShellArray[4] = {oShellSmall,oShellMed,oShellLarge,oShellGiant};
Vector3d latx = Vector3d(lenx,0,0);
Vector3d laty = Vector3d(0,leny,0);
Vector3d latz = Vector3d(0,0,lenz);
Vector3d lat[3] = {latx,laty,latz};
bool flat_div = false; //the divisions will be equal and will divide from z wall to z wall

bool periodic[3] = {false, false, false};
const double dxmin = 0.1;
inline double max(double a, double b) { return (a>b)? a : b; }

int main(int argc, char *argv[]){
  if (argc < 5) {
    printf("usage:  %s Nspheres iterations*N uncertainty_goal filename \n there will be more!\n", argv[0]);
    return 1;
  }

  double maxrad = 0;
  for (int a=5; a<argc; a+=2){
    printf("Checking a = %d which is %s\n", a, argv[a]);
    if (strcmp(argv[a],"outerSphere") == 0) {
      spherical_outer_wall = true;
      periodic[0] = periodic[1] = periodic[2] = false;
      rad = atof(argv[a+1]);
      maxrad = max(maxrad,rad);
      printf("Using outerSphere of %g\n", rad);
    } else if (strcmp(argv[a],"innerSphere") == 0) {
      spherical_inner_wall = true;
      innerRad = atof(argv[a+1]);
    } else if (strcmp(argv[a],"periodxyz") == 0) {
      periodic[0] = true;
      periodic[1] = true;
      periodic[2] = true;
      lenx = atof(argv[a+1]);
      leny = atof(argv[a+1]);
      lenz = atof(argv[a+1]);
      rad = lenx/2;
      maxrad = max(maxrad, rad);
    } else if (strcmp(argv[a],"periodxy") == 0) {
      periodic[0] = true;
      periodic[1] = true;
      lenx = atof(argv[a+1]);
      leny = atof(argv[a+1]);
      rad = lenx/2;
      maxrad = max(maxrad, rad);
    } else if (strcmp(argv[a],"periodx") == 0) {
      periodic[0] = true;
      lenx = atof(argv[a+1]);
      rad = lenx/2;
      maxrad = max(maxrad, rad);
    } else if (strcmp(argv[a],"periody") == 0) {
      periodic[1] = true;
      leny = atof(argv[a+1]);
      maxrad = max(maxrad, leny/2);
    } else if (strcmp(argv[a],"periodz") == 0) {
      periodic[2] = true;
      lenz = atof(argv[a+1]);
      maxrad = max(maxrad, lenz/2);
    } else if (strcmp(argv[a],"wallx") == 0) {
      has_x_wall = true;
      lenx = atof(argv[a+1]);
      periodic[0] = false;
      maxrad = max(maxrad, lenx);
    } else if (strcmp(argv[a],"wally") == 0) {
      has_y_wall = true;
      leny = atof(argv[a+1]);
      periodic[1] = false;
      maxrad = max(maxrad, leny);
    } else if (strcmp(argv[a],"wallz") == 0) {
      has_z_wall = true;
      lenz = atof(argv[a+1]);
      periodic[2] = false;
      maxrad = max(maxrad, lenz);
    } else if (strcmp(argv[a],"flatdiv") == 0) {
      flat_div = true; //otherwise will default to radial divisions
      a -= 1;
    } else {
      printf("Bad argument:  %s\n", argv[a]);
      return 1;
    }
  }
  printf("flatdiv = %s\n", flat_div ? "true" : "false");
  printf("outerSphere = %s\n", spherical_outer_wall ? "true" : "false");
  printf("innerSphere = %s\n", spherical_inner_wall ? "true" : "false");
  latx = Vector3d(lenx,0,0);
  laty = Vector3d(0,leny,0);
  latz = Vector3d(0,0,lenz);
  lat[0] = latx;
  lat[1] = laty;
  lat[2] = latz;


  const char *outfilename = argv[4];
  fflush(stdout);
  const long N = atol(argv[1]);
  const long iterations = long(atol(argv[2])/N*rad*rad*rad/10/10/10);
  const double uncertainty_goal = atof(argv[3]);
  Vector3d *spheres = new Vector3d[N];
  if (uncertainty_goal < 1e-12 || uncertainty_goal > 1.0) {
    printf("Crazy uncertainty goal:  %s\n", argv[1]);
    return 1;
  }
  printf("running with %ld spheres for %ld iterations.\n", N, iterations);

  //////////////////////////////////////////////////////////////////////////////////////////
  // We start with randomly-placed spheres, and then gradually wiggle
  // them around until they are all within bounds and are not
  // overlapping.  We do this by creating an "overlap" value which we
  // constrain to never increase.  Note that this may not work at all
  // for high filling fractions, since we could get stuck in a local
  // minimum.
  for(long i=0; i<N; i++) {
    spheres[i]=10*rad*ran3();
    if (spherical_outer_wall) {
      while (spheres[i].norm() > rad) {
        spheres[i] *= 0.9;
      }
    }
    if (spherical_inner_wall) {
      while (spheres[i].norm() < innerRad) {
        spheres[i] *= rad/innerRad;
      }
    }
  }
  clock_t start = clock();
  long num_to_time = 100000;
  long num_timed = 0;
  long i = 0;
  double scale = .005;

  // Let's move each sphere once, so they'll all start within our
  // periodic cell!
  for (i=0;i<N;i++) spheres[i] = move(spheres[i], scale);

  clock_t starting_initial_state = clock();
  printf("Initial countOverLaps is %g\n", countOverLaps(spheres, N, R));
  while (countOverLaps(spheres, N, R)>0){
    for (int movethis=0;movethis < 100*N; movethis++) {
      if (num_timed++ > num_to_time) {
        clock_t now = clock();
        //printf("took %g seconds per initialising iteration\n",
        //       (now - double(start))/CLOCKS_PER_SEC/num_to_time);
        num_timed = 0;
        start = now;
      }
      Vector3d old =spheres[i%N];
      double oldoverlap = countOneOverLap(spheres, N, i%N, R);
      spheres[i%N]=move(spheres[i%N],scale);
      double newoverlap = countOneOverLap(spheres, N, i%N, R);
      if(newoverlap>oldoverlap){
        spheres[i%N]=old;
      }
      i++;
      if (i%(100*N) == 0) {
        if (i>iterations/4) {
          for(long i=0; i<N; i++) {
            printf("%g\t%g\t%g\n", spheres[i][0],spheres[i][1],spheres[i][2]);
          }
          printf("couldn't find good state\n");
          exit(1);
        }
        char *debugname = new char[10000];
        sprintf(debugname, "%s.debug", outfilename);
        FILE *spheredebug = fopen(debugname, "w");
        for(long i=0; i<N; i++) {
          fprintf(spheredebug, "%g\t%g\t%g\n", spheres[i][0],spheres[i][1],spheres[i][2]);
        }
        fclose(spheredebug);
        printf("numOverLaps=%g (debug file: %s)\n",countOverLaps(spheres,N,R), debugname);
        delete[] debugname;
        fflush(stdout);
      }
    }
  }
  assert(countOverLaps(spheres, N, R) == 0);
  {
    clock_t now = clock();
    //printf("took %g seconds per initialising iteration\n",
    //       (now - double(start))/CLOCKS_PER_SEC/num_to_time);
    printf("\nFound initial state in %g days!\n", (now - double(starting_initial_state))/CLOCKS_PER_SEC/60.0/60.0/24.0);
  }

  long div = uncertainty_goal*uncertainty_goal*iterations;
  if (div < 10) div = 10;
  if (maxrad/div < dxmin) div = int(maxrad/dxmin);
  printf("Using %ld divisions, dx ~ %g\n", div, maxrad/div);
  fflush(stdout);

  double *radius = new double[div+1];
  double *sections = new double [div+1];

  if (flat_div){
    double size = lenz/div;
    for (long s=0; s<div+1; s++){
      sections[s] = size*s - lenz/2.0;
    }
  } else {
    if (spherical_inner_wall){
      double size = rad/div;
      for (long s=0; s<div+1; s++){
        radius[s] = size*s;
      }
    } else {
      const double w = 1.0/(1 + dxmin*div);
      for (long l=0;l<div+1;l++) {
        // make each bin have about the same volume
        radius[l] = w*rad*pow(double(l)/(div), 1.0/3.0) + (1-w)*rad*double(l)/div;
      }
    }
  }

  if (flat_div){
    for (long k=0;k<div+1;k++) printf("section  = %f\n",sections[k]);
    fflush(stdout);
  } else{
    for (long k=0;k<div+1;k++) printf("radius  = %f\n",radius[k]);
    fflush(stdout);
  }

  // Here we use a hokey heuristic to decide on an average move
  // distance, which is proportional to the mean distance between
  // spheres.
  const double mean_spacing = pow(rad*rad*rad/N, 1.0/3);
  if (mean_spacing > 2*R) {
    scale = 2*(mean_spacing - 2*R);
  } else {
    scale = 0.1;
  }
  printf("Using scale of %g\n", scale);
  long count = 0;
  long *shells = new long[div];
  for (long l=0; l<div; l++) shells[l] = 0;

  double *shellsArea = new double [div];
  for (long l=0; l<div; l++) shellsArea[l]=0;
  double *shellsDoubleArea = new double [div];
  for (long l=0; l<div; l++) shellsDoubleArea[l]=0;

  double *density = new double[div];
  double *n0 = new double[div];
  double *nA = new double[div];

  double *SconDensity = new double[div]; double *ScenConDensity = new double[div];
  long *SconShells = new long[div]; long *ScenConShells = new long[div];
  double *MconDensity = new double[div]; double *McenConDensity = new double[div];
  long *MconShells = new long[div]; long *McenConShells = new long[div];
  double *LconDensity = new double[div]; double *LcenConDensity = new double[div];
  long *LconShells = new long[div]; long *LcenConShells = new long[div];
  double *GconDensity = new double[div]; double *GcenConDensity = new double[div];
  long *GconShells = new long[div]; long *GcenConShells = new long[div];


  for(int l=0; l<div; l++){
    SconShells[l]=0; MconShells[l]=0; LconShells[l]=0; GconShells[l]=0;
    ScenConShells[l]=0;McenConShells[l]=0;LcenConShells[l]=0;GcenConShells[l]=0;
  }
  /////////////////////////////////////////////////////////////////////////////

  start = clock();
  num_timed = 0;
  double secs_per_iteration = 0;
  long workingmoves=0;
  int * max_move_counter = new int [N];
  for (int i=0;i<N;i++){max_move_counter[i]=0;}
  int * move_counter = new int [N];
  for (int i=0;i<N;i++){move_counter[i]=0;}
  for (long j=0; j<iterations; j++){
    if (num_timed++ > num_to_time) {
      clock_t now = clock();
      secs_per_iteration = (now - double(start))/CLOCKS_PER_SEC/num_to_time;
      printf("took %g microseconds per iteration\n", 1000000*secs_per_iteration);
      num_timed = 0;
      start = now;
      // after the first timing, just time things once per percent (as
      // often as we print the % complete messages)
      if (iterations/100 > num_to_time) num_to_time = iterations/100;
    }
    // only write out the sphere positions after they've all had a
    // chance to move
    if (workingmoves%N == 0) {
      for (long i=0;i<N;i++) {
        //printf("Sphere at %.1f %.1f %.1f\n", spheres[i][0], spheres[i][1], spheres[i][2]);
        shells[shell(spheres[i], div, radius, sections)]++;
        if (!flat_div){
          for (long k=0; k<div; k++) {
            const double ri = distance(spheres[i],Vector3d(0,0,0));
            // In the following, we compute the shellsArea, which is
            // used to compute n_0, which is one of the fundamental
            // measures of Fundamental Measure Theory.  shellsArea
            // tracks the amount of area of spheres that overlaps with
            // each division of our volume (a.k.a. shell).
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
            // In the following, we will accumulate shellsDoubleArea,
            // which tracks the amount of surface area present in each
            // shell for double-spheres with a radius equal to the
            // *diameter* of a hard sphere.  This is useful in
            // computing the asymmetrically averaged version of the
            // correlation function.
            if (ri < radius[k+1] + 2*R && ri + radius[k+1] > 2*R && ri > radius[k] - 2*R) {
              // There is at least some overlap with shell k! (not so easy)
              double costhetamax, costhetamin;
              if (ri > radius[k] + 2*R) {
                costhetamin = 1;
              } else if (radius[k] + ri < 2*R) {
                costhetamin = 1;
              } else {
                costhetamin = (ri*ri - radius[k]*radius[k] + 2*R*2*R)/(2*ri*2*R);
              }
              if (ri < radius[k+1] - 2*R) {
                costhetamax = -1;
              } else {
                costhetamax = (ri*ri - radius[k+1]*radius[k+1] + 2*R*2*R)/(2*ri*2*R);
              }
              assert(radius[k+1]>radius[k]);
              assert(costhetamin >= costhetamax);
              shellsDoubleArea[k] += 2*M_PI*2*R*2*R*(costhetamin-costhetamax);
            }
          }
        } else {
          // In the following, we compute the shellsArea.  See above
          // for description and discussion.
          for (long k=0; k<div+1; k++){
            double dl = spheres[i][2] - sections[k];
            double dh = spheres[i][2] - sections[k+1];
            if (dl > R) dl = R;
            if (dl < -R) dl = -R;
            if (dh > R) dh = R;
            if (dh < -R) dh = -R;
            shellsArea[k] += 2*M_PI*R*(dl-dh);
          }
          // In the following, we compute the shellsDoubleArea.  See
          // above for description and discussion.
          for (long k=0; k<div+1; k++){
            double dl = spheres[i][2] - sections[k];
            double dh = spheres[i][2] - sections[k+1];
            if (dl > 2*R) dl = 2*R;
            if (dl < -2*R) dl = -2*R;
            if (dh > 2*R) dh = 2*R;
            if (dh < -2*R) dh = -2*R;
            shellsDoubleArea[k] += 2*M_PI*2*R*(dl-dh);
          }
        }
      }
      for(long k=0; k<N; k++){
        for(long n = 0; n<N; n++){
          if (k!=n && touch(spheres[n],spheres[k],oShellArray[3])) {
            GconShells[shell(spheres[k],div,radius,sections)]++;
            GcenConShells[shell((spheres[n]+spheres[k])/2,div,radius,sections)]++;
            if (touch(spheres[n],spheres[k],oShellArray[2])) {
              LconShells[shell(spheres[k],div,radius,sections)]++;
              LcenConShells[shell((spheres[n]+spheres[k])/2,div,radius,sections)]++;
              if (touch(spheres[n],spheres[k],oShellArray[1])) {
                MconShells[shell(spheres[k],div,radius,sections)]++;
                McenConShells[shell((spheres[n]+spheres[k])/2,div,radius,sections)]++;
                if (touch(spheres[n],spheres[k],oShellArray[0])) {
                  SconShells[shell(spheres[k],div,radius,sections)]++;
                  ScenConShells[shell((spheres[n]+spheres[k])/2,div,radius,sections)]++;
                }
              }
            }
          }
        }
      }
    }
    if(j % (iterations/100)==0 && j != 0){
      double secs_to_go = secs_per_iteration*(iterations - j);
      long mins_to_go = secs_to_go / 60;
      long hours_to_go = mins_to_go / 60;
      mins_to_go = mins_to_go % 60;
      if (hours_to_go > 5) {
        printf("%.0f%% complete... (%ld hours to go)\n",j/(iterations*1.0)*100, hours_to_go);
      } else if (mins_to_go < 1) {
        printf("%.0f%% complete... (%.1f seconds to go)\n",j/(iterations*1.0)*100, secs_to_go);
      } else if (hours_to_go < 1) {
        printf("%.0f%% complete... (%ld minutes to go)\n",j/(iterations*1.0)*100, mins_to_go);
      } else if (hours_to_go < 2) {
        printf("%.0f%% complete... (1 hour, %ld minutes to go)\n",j/(iterations*1.0)*100, mins_to_go);
      } else {
        printf("%.0f%% complete... (%ld hours, %ld minutes to go)\n",j/(iterations*1.0)*100, hours_to_go, mins_to_go);
      }
      char *debugname = new char[10000];
      sprintf(debugname, "%s.debug", outfilename);
      FILE *spheredebug = fopen(debugname, "w");
      for(long i=0; i<N; i++) {
        fprintf(spheredebug, "%g\t%g\t%g\n", spheres[i][0],spheres[i][1],spheres[i][2]);
      }
      fclose(spheredebug);
      delete[] debugname;
      fflush(stdout);
    }
    Vector3d temp = move(spheres[j%N],scale);
    count++;
    if(overlap(spheres, temp, N, R, j%N)){
      if (scale > 0.001 && false) {
        scale = scale/sqrt(1.02);
        //printf("Reducing scale to %g\n", scale);
      }
      move_counter[j%N]++;
      if(move_counter[j%N] > max_move_counter[j%N]){
        max_move_counter[j%N] = move_counter[j%N];
      }
      continue;
    }
    move_counter[j%N] = 0;
    spheres[j%N] = temp;
    workingmoves++;
    if (scale < 5 && false) {
      scale = scale*1.02;
      //printf("Increasing scale to %g\n", scale);
    }
  }
  char * counterout = new char[10000];
  sprintf(counterout, "monte-carlo-count-%s-%d.dat", argv[1], int (rad));
  FILE *countout = fopen(counterout,"w");
  for (long i=0;i<N; i++){
    fprintf(countout, "%d\n", max_move_counter[i]);
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  for(long i=0; i<div; i++){
    printf("Number of spheres in division %ld = %ld\n", i+1, shells[i]);
  }
  if (!flat_div){
    for(long i=0; i<div; i++){
      double rmax = radius[i+1];
      double rmin = radius[i];
      density[i]=shells[i]/(((4/3.*M_PI*rmax*rmax*rmax)-(4/3.*M_PI*rmin*rmin*rmin)))/(iterations/double(N));
      n0[i]=shellsArea[i]/(((4/3.*M_PI*rmax*rmax*rmax)-(4/3.*M_PI*rmin*rmin*rmin)))/(iterations/double(N))/(4*M_PI*R*R);
      nA[i]=shellsDoubleArea[i]/(((4/3.*M_PI*rmax*rmax*rmax)-(4/3.*M_PI*rmin*rmin*rmin)))/(iterations/double(N))
        /(4*M_PI*2*R*2*R);
    }
  } else {
    for(long i=0; i<div; i++){
      density[i]=shells[i]/(lenx*leny*lenz/div)/(iterations/double(N));
      n0[i]=shellsArea[i]/(lenx*leny*lenz/div)/(iterations/double(N))/(4*M_PI*R*R);
      nA[i]=shellsDoubleArea[i]/(lenx*leny*lenz/div)/(iterations/double(N))/(4*M_PI*2*R*2*R);
    }
  }

  for(long i=0; i<div; i++){
    SconDensity[i]=((SconShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShellArray[0]*8*oShellArray[0]*oShellArray[0]-4/3.*M_PI*8*R*R*R));
    ScenConDensity[i]=4*M_PI*R*R*((ScenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShellArray[0]*oShellArray[0]*oShellArray[0]-4/3.*M_PI*8*R*R*R));
    MconDensity[i]=((MconShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShellArray[1]*8*oShellArray[1]*oShellArray[1]-4/3.*M_PI*8*R*R*R));
    McenConDensity[i]=4*M_PI*R*R*((McenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShellArray[1]*oShellArray[1]*oShellArray[1]-4/3.*M_PI*8*R*R*R));
    LconDensity[i]=((LconShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShellArray[2]*8*oShellArray[2]*oShellArray[2]-4/3.*M_PI*8*R*R*R));
    LcenConDensity[i]=4*M_PI*R*R*((LcenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShellArray[2]*oShellArray[2]*oShellArray[2]-4/3.*M_PI*8*R*R*R));
    GconDensity[i]=((GconShells[i]+0.0)/shells[i])/((4/3.*M_PI*oShellArray[3]*8*oShellArray[3]*oShellArray[3]-4/3.*M_PI*8*R*R*R));
    GcenConDensity[i]=4*M_PI*R*R*((GcenConShells[i]+0.0)/shellsArea[i])/((4/3.*M_PI*8*oShellArray[3]*oShellArray[3]*oShellArray[3]-4/3.*M_PI*8*R*R*R));
  }
  for(long i=0; i<div; i++){
    printf("Number of contacts in division %ld = %ld %ld %ld %ld\n", i+1,
           SconShells[i], MconShells[i], LconShells[i], GconShells[i]);
    printf("Number of contacts (center) in division %ld = %ld %ld %ld %ld\n", i+1,
           ScenConShells[i], McenConShells[i], LcenConShells[i], GcenConShells[i]);
  }

  //FILE *out = fopen((const char *)outfilename,"w");
  FILE *out = fopen((const char *)outfilename,"w");
  if (out == NULL) {
    printf("Error creating file %s\n", outfilename);
    return 1;
  }
  if (flat_div){
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 0.5*(sections[0]+sections[1]), density[0],
            SconDensity[0], ScenConDensity[0], MconDensity[0], McenConDensity[0],
            LconDensity[0], LcenConDensity[0], GconDensity[0], GcenConDensity[0], n0[0], nA[0]);
  } else if (spherical_inner_wall) {
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", radius[0], 0.0,
            SconDensity[0], ScenConDensity[0], MconDensity[0], McenConDensity[0],
            LconDensity[0], LcenConDensity[0], GconDensity[0], GcenConDensity[0], n0[0], nA[0]);
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", 0.5*(radius[0]+radius[1]), density[0],
            SconDensity[0], ScenConDensity[0], MconDensity[0], McenConDensity[0],
            LconDensity[0], LcenConDensity[0], GconDensity[0], GcenConDensity[0], n0[0], nA[0]);
  } else {
    fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n" , 0.0, density[0],
            SconDensity[0], ScenConDensity[0], MconDensity[0], McenConDensity[0],
            LconDensity[0], LcenConDensity[0], GconDensity[0], GcenConDensity[0], n0[0], nA[0]);
  }

  long divtoprint = div;
  if (!spherical_outer_wall) divtoprint = div - 1;
  if (!flat_div) {
    for(long i=1; i<divtoprint; i++){
      fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
              0.5*(radius[i]+radius[i+1]), density[i],
              SconDensity[i], ScenConDensity[i], MconDensity[i], McenConDensity[i],
              LconDensity[i], LcenConDensity[i], GconDensity[i], GcenConDensity[i], n0[i], nA[i]);
    }
  } else {
    for(long i=1; i<div; i++){
      fprintf(out, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
              0.5*(sections[i]+sections[i+1]), density[i],
              SconDensity[i], ScenConDensity[i], MconDensity[i], McenConDensity[i],
              LconDensity[i], LcenConDensity[i], GconDensity[i], GcenConDensity[i], n0[i], nA[i]);
    }
  }

  printf("Total number of attempted moves = %ld\n",count);
  printf("Total number of successful moves = %ld\n",workingmoves);
  printf("Acceptance rate = %g\n", workingmoves/double(count));
  //delete[] shells;
  fflush(stdout);
  //delete[] density;
  fflush(stdout);
  delete[] SconShells;
  fflush(stdout);
  delete[] SconDensity;
  fflush(stdout);
  delete[] ScenConDensity; delete[] ScenConShells;
  delete[] MconShells; delete[] MconDensity; delete[] McenConDensity; delete[] McenConShells;
  delete[] LconShells; delete[] LconDensity; delete[] LcenConDensity; delete[] LcenConShells;
  delete[] GconShells; delete[] GconDensity; delete[] GcenConDensity; delete[] GcenConShells;
  fflush(stdout);
  delete[] spheres;
  delete[] max_move_counter;
  delete[] move_counter;
  fflush(stdout);
  fclose(out);
  fclose(countout);
}

double countOneOverLap(Vector3d *spheres, long n, long j, double R){
  double num = 0;
  for(long i = 0; i < n; i++){
    if(i != j && distance(spheres[i],spheres[j])<2*R){
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
      num += innerRad - distance(spheres[j],Vector3d(0,0,0));
    }
  }
  if (has_x_wall || periodic[0]){
    if (spheres[j][0] > lenx/2 ){
      num += spheres[j][0]-(lenx/2);
    } else if (spheres[j][0] < -lenx/2){
      num -= spheres[j][0] + (lenx/2);
    }
  }
  if (has_y_wall || periodic[1]){
    if (spheres[j][1] > leny/2 ){
      num += spheres[j][1]-(leny/2);
    } else if (spheres[j][1] < -leny/2){
      num -= spheres[j][1] + (leny/2);
    }
  }
  if (has_z_wall || periodic[2]){
    if (spheres[j][2] > lenz/2 ){
      num += spheres[j][2]-(lenz/2);
    } else if (spheres[j][2] < -lenz/2){
      num -= spheres[j][2] + (lenz/2);
    }
  }
  for(long i = 0; i < n; i++){
    if (i != j) {
      for (long k=0; k<3; k++){
	if (periodic[k]){
	  if (distance(spheres[j],spheres[i]+lat[k]) < 2*R){
	    num += 2*R - distance(spheres[j],spheres[i]+lat[k]);
	  } else if (distance(spheres[j],spheres[i]-lat[k]) < 2*R) {
	    num += 2*R - distance(spheres[j],spheres[i]-lat[k]);
	  }
	}
	for (long m=k+1; m<3; m++){
	  if (periodic[m] && periodic[k]){
	    if (distance(spheres[j],spheres[i]+lat[k]+lat[m]) < 2*R){
	      num += 2*R - distance(spheres[j],spheres[i]+lat[k]+lat[m]);
	    }
	    if (distance(spheres[j],spheres[i]-lat[k]-lat[m]) < 2*R){
	      num += 2*R - distance(spheres[j],spheres[i]-lat[k]-lat[m]);
	    }
	    if (distance(spheres[j],spheres[i]+lat[k]-lat[m]) < 2*R){
	      num += 2*R - distance(spheres[j],spheres[i]+lat[k]-lat[m]);
	    }
            if (distance(spheres[j],spheres[i]-lat[k]+lat[m]) < 2*R){
              num += 2*R - distance(spheres[j],spheres[i]-lat[k]+lat[m]);
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



double countOverLaps(Vector3d *spheres, long n, double R){
  double num = 0;
  for(long j = 0; j<n; j++){
    num += countOneOverLap(spheres, n, j, R);
  }
  return num;
}

bool overlap(Vector3d *spheres, Vector3d v, long n, double R, long s){
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
  bool amonborder[3] = {
    fabs(v[0]) + 2*R >= lenx/2,
    fabs(v[1]) + 2*R >= leny/2,
    fabs(v[2]) + 2*R >= lenz/2
  };
  for(long i = 0; i < n; i++){
    if (i!=s){
      if(distance(spheres[i],v)<2*R){
        return true;
      }
    }
  }
  for (long k=0; k<3; k++) {
    if (periodic[k] && amonborder[k]) {
      for(long i = 0; i < n; i++) {
        if (i!=s){
          if (distance(v,spheres[i]+lat[k]) < 2*R) return true;
          if (distance(v,spheres[i]-lat[k]) < 2*R) return true;
        }
      }
      for (long m=k+1; m<3; m++){
        if (periodic[m] && amonborder[m]){
          for(long i = 0; i < n; i++) {
            if (i!=s){
              if (distance(v,spheres[i]+lat[k]+lat[m]) < 2*R) return true;
              if (distance(v,spheres[i]-lat[k]-lat[m]) < 2*R) return true;
              if (distance(v,spheres[i]+lat[k]-lat[m]) < 2*R) return true;
              if (distance(v,spheres[i]-lat[k]+lat[m]) < 2*R) return true;
            }
          }
        }
      }
    }
  }
  if (periodic[0] && periodic[1] && periodic[2]
      && amonborder[0] && amonborder[1] && amonborder[2]){
    for(long i = 0; i < n; i++) {
      if (i!=s){
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
  if (periodic[0] || has_x_wall){
    while (newv[0] > lenx/2){
      newv[0] -= lenx;
    }
    while (newv[0] < -lenx/2){
      newv[0] += lenx;
    }
  }
  if (periodic[1] || has_y_wall){
    while (newv[1] > leny/2){
      newv[1] -= leny;
    }
    while (newv[1] < -leny/2){
      newv[1] += leny;
    }
  }
  if (periodic[2] || has_z_wall){
    while (newv[2] > lenz/2){
      newv[2] -= lenz;
    }
    while (newv[2] < -lenz/2){
      newv[2] += lenz;
    }
  }
  //printf("Moved to %.1f %.1f %.1f by scale %g\n", newv[0], newv[1], newv[2], scale);
  return newv;
}


bool touch(Vector3d w, Vector3d v, double oShell){
  const double dvw = distance(v,w);
  if (dvw < 2*oShell) return true;
  // The following is a hack to avoid doing many distance calculations
  // in cases where we can be certain that periodic boundaries
  // couldn't be allowing these two spheres to touch.  It makes a
  // shocking difference in the overall speed, when we have periodic
  // boundary conditions in all three directions!
  if (dvw < lenx - 2*oShell &&
      dvw < leny - 2*oShell &&
      dvw < lenz - 2*oShell) return false;

  // Now we check for all the possible ways the two spheres could
  // touch across the cell in periodic directions...
  //return false;
  for (long k=0; k<3; k++){
    if (periodic[k]){
      if (distance(v,w+lat[k]) < 2*oShell) return true;
      if (distance(v,w-lat[k]) < 2*oShell) return true;
      for (long m=k+1; m<3; m++){
        if (periodic[m]){
          if (distance(v,w+lat[k]+lat[m]) < 2*oShell) return true;
          if (distance(v,w-lat[k]-lat[m]) < 2*oShell) return true;
          if (distance(v,w+lat[k]-lat[m]) < 2*oShell) return true;
          if (distance(v,w-lat[k]+lat[m]) < 2*oShell) return true;
        }
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



long shell(Vector3d v, long div, double *radius, double *sections){
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
  const long unsigned int x =0;
  static MTRand my_mtrand(x); // always use the same random number generator (for debugging)!
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
