#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "utilities.h"
using std::vector;

inline double countOverLaps(Vector3d *spheres, long n, double R);
inline double countOneOverLap(Vector3d *spheres, long n, long j, double R);
inline bool overlap(Vector3d *spheres, Vector3d v, long n, double R, long s);
inline double distXY(Vector3d a, Vector3d b);
inline double distXYZ(Vector3d a, Vector3d b);
inline Vector3d periodicDiff(Vector3d a, Vector3d b);
static void took(const char *name);

const double dz = 0.1;
const double dx = 0.1;

const double path_dr = 0.01;
const double path_dx = 0.1;
const double path_dtheta = 2.0/3.0*M_PI/42.0;

// resolution info for the a1 histogram / integral
const double a1_dr = 0.01;
const double a1_dz = 0.01;
const double a1_rmax = 5.0; // max value to store for delta x

bool path = false; // only look at where |z1-z0| = 2 and |z1-z0| = 4
                   // which are the distances where we store path info
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
Vector3d latx = Vector3d(lenx,0,0);
Vector3d laty = Vector3d(0,leny,0);
Vector3d latz = Vector3d(0,0,lenz);
Vector3d lat[3] = {latx,laty,latz};
bool flat_div = true; // the divisions will be equal and will divide from z wall to z wall

bool periodic[3] = {false, false, false};
const double dxmin = 0.1;
inline double max(double a, double b) { return (a>b)? a : b; }
inline double min(double a, double b) { return (a<b)? a : b; }

int main(int argc, char *argv[]){
  if (argc < 5) {
    printf("usage:  %s Nspheres iterations*N uncertainty_goal filename a1_filename\n there will be more!\n", argv[0]);
    return 1;
  }
  double maxrad = 0;
  for (int a=6; a<argc; a+=2){
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
    } else if (strcmp(argv[a],"path") == 0) {
      path = true;
      a -= 1;
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

  const int zbins = path ? 2 : lenx/2.0/dz;
  const int xbins = lenx/2.0/dx;
  const int z2bins = lenx/2.0/dz;

  const int a1_zbins = lenx/2.0/a1_dz;
  const int a1_rbins = path ? 2 : (a1_rmax - 2.0)/a1_dr;

  const int path_zbins = (lenx/2.0 - 4*R)/path_dr;
  const int path_thetabins = 2.0/3.0*M_PI/path_dtheta;
  const int path_xbins = (lenx/2.0 - R)/path_dr;

  const int path2_zbins = (lenx/2.0 - 6*R)/path_dr;
  const int path2_thetabins = M_PI/path_dtheta;
  const int path2_xbins = (lenx/2.0)/path_dr;

  const char *outfilename = argv[4];
  const char *da_dz_outfilename = argv[5];
  char *pathfilename = new char[1024];
  fflush(stdout);
  const long N = atol(argv[1]);
  const double density = N/lenx/leny/lenz;
  const long iterations = long(atol(argv[2])/N*rad*rad*rad/10/10/10);
  const double uncertainty_goal = atof(argv[3]);
  Vector3d *spheres = new Vector3d[N];


  long *histogram = new long[zbins*xbins*z2bins]();
  long *da_dz_histogram = new long[a1_rbins*a1_zbins]();
  const int path_bins = path_xbins + path_thetabins + path_zbins;
  const int path2_bins = path2_xbins + path2_thetabins + path2_zbins;
  long *path_histogram = new long[path_bins]();
  long *path2_histogram = new long[path2_bins]();
  long numinhistogram = 0;
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
  long num_to_time = 100*N;
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
      Vector3d old = spheres[i%N];
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
        if (!spheredebug) {
          printf("oops opening %s\n", debugname);
          exit(1);
        }
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

  /////////////////////////////////////////////////////////////////////////////
  int hours_now = 1;
  num_to_time = 1000;
  start = clock();
  num_timed = 0;
  double secs_per_iteration = 0;
  long workingmoves=0;
  int * max_move_counter = new int [N];
  for (int i=0;i<N;i++){max_move_counter[i]=0;}
  int * move_counter = new int [N];
  for (int i=0;i<N;i++){move_counter[i]=0;}

  clock_t output_period = CLOCKS_PER_SEC*60; // start at outputting every minute
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*60; // top out at one hour interval
  clock_t last_output = clock(); // when we last output data
  took("Setting up");
  // Begin main program loop
  //////////////////////////////////////////////////////////////////////////////
  for (long j=0; j<iterations; j++){
    // if (j%(N*1000) == 0) {
    //   took("One thousand iterations");
    //   printf("count: %li\n", count);
    // }
	  num_timed = num_timed + 1;
    if (num_timed > num_to_time || j==(iterations-1)) {
      num_timed = 0;
      ///////////////////////////////////////////start of print.dat
      const clock_t now = clock();
      secs_per_iteration = (now - double(start))/CLOCKS_PER_SEC/num_to_time;
      if (secs_per_iteration*num_to_time < 1) {
        printf("took %g microseconds per iteration\n", 1000000*secs_per_iteration);
        num_to_time *= 2;
      } else {
        // Set the number of iterations to time to a minute, so we
        // won't check *too* many times.
        num_to_time = long(60/secs_per_iteration);
      }
      start = now;
      if ((now > last_output + output_period) || j==(iterations-1)) {
        last_output = now;
        if (output_period < max_output_period/2) {
          output_period *= 2;
        } else if (output_period < max_output_period) {
          output_period = max_output_period;
        }
        {
          double secs_done = double(now)/CLOCKS_PER_SEC;
          long mins_done = secs_done / 60;
          long hours_done = mins_done / 60;
          mins_done = mins_done % 60;
          if (hours_done > 50) {
            printf("Saved data after %ld hours\n", hours_done);
          } else if (mins_done < 1) {
            printf("Saved data after %.1f seconds\n", secs_done);
          } else if (hours_done < 1) {
            printf("Saved data after %ld minutes\n", mins_done);
          } else if (hours_done < 2) {
            printf("Saved data after 1 hour, %ld minutes\n", mins_done);
          } else {
            printf("Saved data after %ld hours, %ld minutes\n", hours_done, mins_done);
          }
        }
        if (flat_div){
          // save the da_dz data
          char *da_dz_filename = new char[1024];
          for (int i=0; i<a1_rbins; i++) {
            const double a1_r12 = path ? 2*R*(i+1) + 0.5*a1_dr : 2*R + (i+0.5)*a1_dr;
            sprintf(da_dz_filename, "%s-%5.3f.dat", da_dz_outfilename, a1_r12);
            FILE *da_dz_out = fopen((const char *)da_dz_filename, "w");
            for (int k=0; k<a1_zbins; k++) {
              const double a1_z1_max = (k+1)*a1_dr;
              const double a1_z1_min = k*a1_dr;
              const double slice_volume = 4.0/3.0*M_PI*(a1_z1_max*a1_z1_max*a1_z1_max -
                                                       a1_z1_min*a1_z1_min*a1_z1_min);
              const double da_dz = double(da_dz_histogram[i*a1_zbins + k])/
                slice_volume/a1_dr/double(count);
              const double coord = (k + 0.5)*a1_dz;
              fprintf(da_dz_out, "%g\t%g\n", coord, da_dz);
            }
            fclose(da_dz_out);
          }
          delete[] da_dz_filename;
          // // save the path data
          sprintf(pathfilename, "%s-path.dat", outfilename);
          FILE *path_out = fopen((const char *)pathfilename, "w");
          if (path_out == NULL) {
            printf("Error creating file %s\n", pathfilename);
            return 1;
          }
          fprintf(path_out, "# Working moves: %li, total moves: %li\n", workingmoves, count);
          fprintf(path_out, "# s       g3         z         x         vol2        histogram\n");
          const double path_vol0 = lenx*leny*lenz;
          const double path_r1max = 2.0*R + path_dr;
          const double path_r1min = 2.0*R;
          const double path_vol1 = 4.0/3.0*M_PI*(path_r1max*path_r1max*path_r1max -
                                                 path_r1min*path_r1min*path_r1min);
          double s = -path_dr/2.0;
          for(int i=0; i<path_zbins; i++) {
            const double probability = double(path_histogram[i])
              /double(numinhistogram);
            const double path_vol2 = M_PI*path_dx*path_dx*path_dr;
            const double n3 = probability/path_vol0/path_vol1/path_vol2;
            const double g3 = n3/density/density/density;
            s += path_dr;
            const double z2 = (path_zbins - i)*path_dr + 4*R;
            const double x2 = path_dx/2.0;
            fprintf(path_out, "%6.3f    %7.4f    %6.3f    %6.3f    %8.5f    %li # horizontal leg\n",
                    s,g3,z2,x2,path_vol2,path_histogram[i]);
          }
          // finding the difference between the last point coming over and the first
          // point in the theta band
          const double z0 = path_dr + 4*R;
          const double x0 = path_dx/2.0;
          const double z1 = (2.0*R + path_dr/2.0) + (2.0*R + path_dr/2.0)*cos(path_dtheta/2);
          const double x1 = (2.0*R + path_dr/2.0)*sin(path_dtheta/2);
          s += sqrt((x1-x0)*(x1-x0) + (z1-z0)*(z1-z0));
          ////
          s -= (2.0 + path_dr/2)*path_dtheta;
          for(int i=path_zbins; i<path_zbins+path_thetabins; i++) {
            const double probability = double(path_histogram[i])
              /double(numinhistogram)/2.0;
            const double rmax = 2.0*R + path_dr;
            const double rmin = 2.0*R;
            const double theta = (i - path_zbins + 0.5)*path_dtheta;
            const double theta_min = theta - path_dtheta/2.0;
            const double theta_max = theta + path_dtheta/2.0;
            const double path_vol2 = 2.0/3.0*M_PI*(rmax*rmax*rmax - rmin*rmin*rmin)*
              (cos(theta_min) - cos(theta_max));
            const double n3 = probability/path_vol0/path_vol1/path_vol2;
            const double g3 = n3/density/density/density;
            s += (2.0 + path_dr/2)*path_dtheta;

            const double z2 = (2*R + path_dr/2.0) + (2*R + path_dr/2.0)*cos(theta);
            const double x2 = (2*R + path_dr/2.0)*sin(theta);
            fprintf(path_out, "%6.3f    %7.4f    %6.3f    %6.3f    %8.5f    %li # theta leg\n",
                    s,g3,z2,x2,path_vol2,path_histogram[i]);
          }
          // finding the difference between the last point in the theta band and the
          // first point in the vertical band
          const double theta_m = (path_thetabins - 0.5)*path_dtheta;
          const double z2 = (2.0 + path_dr/2.0)*(1.0 + cos(theta_m));
          const double x2 = (2*R + path_dr/2.0)*sin(theta_m);
          const double z3 = (2*R + path_dr/2.0)/2.0;
          const double x3 = sqrt(3)*R + path_dr/2.0;
          s += sqrt((z3-z2)*(z3-z2) + (x3-x2)*(x3-x2));
          ////
          s -= path_dr;
          for(int i=path_zbins+path_thetabins; i<path_zbins+path_thetabins+path_xbins; i++){
            const double probability = double(path_histogram[i])
              /double(numinhistogram)/2.0;
            const double x2 = (i-path_zbins-path_thetabins+0.5)*path_dr + sqrt(3)*R;
            const double rmax = x2 + path_dr/2.0;
            const double rmin = x2 - path_dr/2.0;
            const double path_vol2 = M_PI*(rmax*rmax - rmin*rmin)*path_dr;
            const double n3 = probability/path_vol0/path_vol1/path_vol2;
            const double g3 = n3/density/density/density;
            s += path_dr;
            const double z2 = (2*R + path_dr/2.0)/2.0;
            fprintf(path_out, "%6.3f    %7.4f    %6.3f    %6.3f    %8.5f    %li # vertical leg\n",
                    s,g3,z2,x2,path_vol2,path_histogram[i]);
          }
          fclose(path_out);
          // // ------------------

          // save the path2 data ----------------------------------------------
          sprintf(pathfilename, "%s-path2.dat", outfilename);
          FILE *path2_out = fopen((const char *)pathfilename, "w");
          if (path2_out == NULL) {
            printf("Error creating file %s\n", pathfilename);
            return 1;
          }
          fprintf(path2_out, "# Working moves: %li, total moves: %li\n", workingmoves, count);
          fprintf(path_out, "# s       g3         z         x         vol2        histogram\n");
          const double path2_vol0 = lenx*leny*lenz;
          const double path2_r1max = 4.0*R + path_dr;
          const double path2_r1min = 4.0*R;
          const double path2_vol1 = 4.0/3.0*M_PI*(path2_r1max*path2_r1max*path2_r1max -
                                                 path2_r1min*path2_r1min*path2_r1min);
          s = -path_dr/2.0;
          for(int i=0; i<path2_zbins; i++) {
            const double probability = double(path2_histogram[i])
              /double(numinhistogram);
            const double path2_vol2 = M_PI*path_dx*path_dx*path_dr;
            const double n3 = probability/path2_vol0/path2_vol1/path2_vol2;
            const double g3 = n3/density/density/density;
            s += path_dr;
            const double z2 = (path2_zbins - i)*path_dr + 6*R;
            const double x2 = path_dx/2.0;
            fprintf(path_out, "%6.3f    %7.4f    %6.3f    %6.3f    %8.5f    %li # horizontal leg\n",
                    s, g3, z2, x2, path2_vol2, path2_histogram[i]);
          }
          // finding the difference between the last point coming over and the first
          // point in the theta band
          const double z4 = path_dr + 6*R;
          const double x4 = path_dx/2.0;
          const double z5 = (4.0*R + path_dr/2.0) + (2.0*R + path_dr/2.0)*cos(path_dtheta/2.0);
          const double x5 = (2.0*R + path_dr/2.0)*sin(path_dtheta/2.0);
          s += sqrt((x5-x4)*(x5-x4) + (z5-z4)*(z5-z4));
          ////
          s -= (2.0 + path_dr/2.0)*path_dtheta;
          for(int i=path2_zbins; i<path2_zbins+path2_thetabins; i++) {
            const double probability = double(path2_histogram[i])
              /double(numinhistogram)/2.0;
            const double rmax = 2.0*R + path_dr;
            const double rmin = 2.0*R;
            const double theta = (i - path2_zbins + 0.5)*path_dtheta;
            const double theta_min = theta - path_dtheta/2.0;
            const double theta_max = theta + path_dtheta/2.0;
            const double path2_vol2 = 2.0/3.0*M_PI*(rmax*rmax*rmax - rmin*rmin*rmin)*
              (cos(theta_min) - cos(theta_max));
            const double n3 = probability/path2_vol0/path2_vol1/path2_vol2;
            const double g3 = n3/density/density/density;
            s += (2.0 + path_dr/2.0)*path_dtheta;

            const double z2 = (4*R + path_dr/2.0) + (2*R + path_dr/2.0)*cos(theta);
            const double x2 = (2*R + path_dr/2.0)*sin(theta);
            fprintf(path_out, "%6.3f    %7.4f    %6.3f    %6.3f    %8.5f    %li # theta leg\n",
                    s, g3, z2, x2, path2_vol2, path2_histogram[i]);
          }
          // finding the difference between the last point in the theta band and the
          // first point in the vertical band
          const double theta_max = (path2_thetabins - 0.5)*path_dtheta;
          const double z6 = (4*R + path_dr/2.0) + (2*R + path_dr/2.0)*cos(theta_max);
          const double x6 = (2*R + path_dr/2.0)*sin(theta_max);
          const double z7 = (2*R + path_dr);
          const double x7 = -path_dr/2.0; // minus sign so the distance is as though
          // we circle around z1 and go through, but we don't actually have do do that
          // because of symmetry
          s += sqrt((z7-z6)*(z7-z6) + (x7-x6)*(x7-x6));
          ////
          s -= path_dr;
          for(int i=path2_zbins+path2_thetabins; i<path2_zbins+path2_thetabins+path2_xbins; i++){
            const double probability = double(path2_histogram[i])
              /double(numinhistogram)/2.0;
            const double x2 = (i-path2_zbins-path2_thetabins+0.5)*path_dr;
            const double rmax = x2 + path_dr/2.0;
            const double rmin = x2 - path_dr/2.0;
            const double path2_vol2 = M_PI*(rmax*rmax - rmin*rmin)*path_dr;
            const double n3 = probability/path2_vol0/path2_vol1/path2_vol2;
            const double g3 = n3/density/density/density;
            s += path_dr;
            const double z2 = (2*R + path_dr);
            fprintf(path_out, "%6.3f    %7.4f    %6.3f    %6.3f    %8.5f    %li # vertical leg\n",
                    s, g3, z2, x2, path2_vol2, path2_histogram[i]);
          }
          fclose(path2_out);
          // --------------------------------------------------------------------

          // save the triplet distribution data
          const double volume0 = lenx*leny*lenz;
          for (int z1_i=0; z1_i<zbins; z1_i++) {
            const double z1 = path ? 2*R*(z1_i+1) + dz/2.0 : 2*R + (z1_i + 0.5)*dz;
            char *finalfilename = new char[1024];
            sprintf(finalfilename, "%s-%05.2f.dat", outfilename, z1);
            FILE *out = fopen((const char *)finalfilename, "w");
            if (out == NULL) {
              printf("Error creating file %s\n", finalfilename);
              return 1;
            }
            delete[] finalfilename;
            const double r1min = z1-dz/2.0;
            const double r1max = z1+dz/2.0;
            const double volume1 = 4.0/3.0*M_PI*(r1max*r1max*r1max - r1min*r1min*r1min);
            for (double x2=dx/2.0; x2<lenx/2.0; x2+=dx) {
              for (double z2=dz/2.0; z2<lenx/2.0; z2+=dz) {
                const double x2min = x2-dx/2.0;
                const double x2max = x2+dx/2.0;
                const double volume2 = M_PI*(x2max*x2max - x2min*x2min)*dz;
                const double probability =
                  double(histogram[z1_i*xbins*z2bins + int(x2/dx)*z2bins +
                                   int(z2/dz)])/double(numinhistogram);
                const double n3 = probability/volume0/volume1/volume2;
                const double g3 = n3/density/density/density;
                fprintf(out, "%g\t", g3);
              }
              fprintf(out, "\n");
            }
          fclose(out);
          }
        }
        fflush(stdout);
      }
      ///////////////////////////////////////////end of print.dat
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

    // only write out the sphere positions after they've all had a
    // chance to move
    // positions are converted from coordinates to matrix indices,
    // then histogram counts spheres at the
    // appropriate locations
    if (count%N == 0) {
      numinhistogram++;
      for (int l=0; l<N; l++) {
        for (int i=0; i<N; i++) {
          if (i != l) {
            const Vector3d r01 = periodicDiff(spheres[i], spheres[l]);
            const double z1 = r01.norm();
            int which_path = -1; // -1 if no path, 0 if touching path, 1 if path between
            if (z1 < 2*R + path_dr) which_path = 0;
            else if (z1 > 4*R && z1 < 4*R + path_dr) which_path = 1;
            if(!path || which_path > -1) {
              const Vector3d zhat = r01.normalized();
              const int z1_i = path ? which_path : int(z1/dz);
              const int a1_z1_i = int(z1/a1_dz);
              //printf("Sphere at %.1f %.1f %.1f\n", spheres[i][0], spheres[i][1], spheres[i][2]);
              for (int k=0; k<N; k++) {
                if (k != i && k != l) { // don't look at a sphere in relation to itself
                  const Vector3d r02 = periodicDiff(spheres[k], spheres[l]);
                  const double z2 = r02.dot(zhat);
                  const int z2_i = floor(z2/dz);
                  const double x2 = (r02 - z2*zhat).norm();
                  const int x2_i = int(x2/dx);

                  const Vector3d r12_v = periodicDiff(spheres[k], spheres[i]);
                  const double r12 = r12_v.norm();
                  const int a1_r12_i = path ? which_path : int((r12-2.0)/a1_dr);

                  if (z1_i < zbins && x2_i < xbins && z2_i >= 0
                      && z2_i < z2bins && r02.norm() < lenx/2.0)
                    histogram[z1_i*xbins*z2bins + x2_i*z2bins + z2_i] ++;

                  if (a1_r12_i < a1_rbins && a1_z1_i < a1_zbins)
                    da_dz_histogram[a1_r12_i*a1_zbins + a1_z1_i] ++;

                  const double z_cent = z1/2.0;
                  const double theta = acos(r12_v.dot(zhat)/r12);
                  if (which_path == 0 && z2 > 0) {
                    if (z2 < lenx/2.0 && x2 < path_dx && z2 > z1 + 2*R) {
                      const int index = path_zbins - ceil((z2-z1-2*R)/path_dr);
                      if (index < 0 || index >= path_zbins)
                        fprintf(stderr,"Index out of bounds: %i, z1: %.2f, z2: %.2f, x2: %.2f\n",index,z1,z2,x2);
                      path_histogram[index] ++;
                    } if (r12 < 2*R + path_dr && theta < 2.0/3.0*M_PI) {
                      const int index = path_zbins + theta/path_dtheta;
                      if (index < path_zbins || index >= path_zbins + path_thetabins)
                        fprintf(stderr,"Index out of bounds: %i, z1: %.2f, z2: %.2f, x2: %.2f\n",index,z1,z2,x2);
                      path_histogram[index] ++;
                    } if (fabs(z2-z_cent) < path_dr && x2 < lenx/2.0) {
                      const double xmin = sqrt(3)*R;
                      const int index = path_zbins + path_thetabins + int((x2 - xmin)/path_dr);
                      if (index < path_zbins+path_thetabins || index >= path_bins)
                        fprintf(stderr,"Index out of bounds: %i, z1: %.2f, z2: %.2f, x2: %.2f\n",index,z1,z2,x2);
                      path_histogram[index] ++;
                    }
                  }
                  else if (which_path == 1 && z2 > 0) {
                    if (z2 < lenx/2.0 && x2 < path_dx && z2 > z1 + 2*R) {
                      const int index = path2_zbins - ceil((z2-z1-2*R)/path_dr);
                      if (index < 0 || index >= path_zbins)
                        fprintf(stderr,"Index2 out of bounds: %i, z1: %.2f, z2: %.2f, x2: %.2f\n",index,z1,z2,x2);
                      path2_histogram[index] ++;
                    } if (r12 < 2*R + path_dr) {
                      const int index = path2_zbins + theta/path_dtheta;
                      if (index < path2_zbins || index >= path2_zbins + path2_thetabins)
                        fprintf(stderr,"Index2 out of bounds: %i, z1: %.2f, z2: %.2f, x2: %.2f\n",index,z1,z2,x2);
                      path2_histogram[index] ++;
                    } if (fabs(z2-z_cent) < path_dr && x2 < lenx/2.0) {
                      const int index = path2_zbins + path2_thetabins + int(x2/path_dr);
                      if (index < path2_zbins+path2_thetabins || index >= path2_bins)
                        fprintf(stderr,"Index2 out of bounds: %i, z1: %.2f, z2: %.2f, x2: %.2f\n",index,z1,z2,x2);
                      path2_histogram[index] ++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
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
    }

  //////////////////////////////////////////////////////////////////////////////
  // End main program loop
  if (scale < 5 && false) {
    scale = scale*1.02;
    //printf("Increasing scale to %g\n", scale);
  }
  char * counterout = new char[10000];
  sprintf(counterout, "monte-carlo-count-%s-%d.dat", argv[1], int (rad));
  FILE *countout = fopen(counterout,"w");
  if (!countout) {
    printf("error opening %s\n", counterout);
    exit(1);
  }
  delete[] counterout;
  for (long i=0;i<N; i++){
    fprintf(countout, "%d\n", max_move_counter[i]);
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  for(long i=0; i<div; i++){
    printf("Number of spheres in division %ld = %ld\n", i+1, shells[i]);
  }

  printf("Total number of attempted moves = %ld\n",count);
  printf("Total number of successful moves = %ld\n",workingmoves);
  printf("Acceptance rate = %g\n", workingmoves/double(count));
  fflush(stdout);
  delete[] spheres;
  delete[] max_move_counter;
  delete[] move_counter;
  delete[] histogram;
  delete[] da_dz_histogram;
  delete[] path_histogram;
  fflush(stdout);
  fclose(countout);
}

inline double countOneOverLap(Vector3d *spheres, long n, long j, double R){
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

inline double countOverLaps(Vector3d *spheres, long n, double R){
  double num = 0;
  for(long j = 0; j<n; j++){
    num += countOneOverLap(spheres, n, j, R);
  }
  return num;
}

inline bool overlap(Vector3d *spheres, Vector3d v, long n, double R, long s){
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

inline Vector3d fixPeriodic(Vector3d newv){
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

inline Vector3d move(Vector3d v, double scale){
  Vector3d newv = v+scale*ran3();
  return fixPeriodic(newv);
}

inline double distXY(Vector3d a, Vector3d b){
  const double xdist = periodic[0] ?
    min(fabs(b.x()-a.x()), (lenx-fabs(b.x()-a.x()))) : b.x()-a.x();
  const double ydist = periodic[1] ?
    min(fabs(b.y()-a.y()), (leny-fabs(b.y()-a.y()))) : b.y()-a.y();
  return sqrt(xdist*xdist + ydist*ydist);
}
inline double distXYZ(Vector3d a, Vector3d b){
  const double xdist = periodic[0] ?
    min(fabs(b.x()-a.x()), (lenx-fabs(b.x()-a.x()))) : b.x()-a.x();
  const double ydist = periodic[1] ?
    min(fabs(b.y()-a.y()), (leny-fabs(b.y()-a.y()))) : b.y()-a.y();
  const double zdist = periodic[2] ?
    min(fabs(b.z()-a.z()), (lenz-fabs(b.z()-a.z()))) : b.z()-a.z();
  return sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
}

inline Vector3d periodicDiff(Vector3d b, Vector3d a){
  Vector3d v = b - a;
  if (periodic[0]) {
    while (v.x() > lenx/2.0)
      v[0] -= lenx;
    while (v.x() < -lenx/2.0)
      v[0] += lenx;
  }
  if (periodic[1]) {
    while (v.y() > leny/2.0)
      v[1] -= leny;
    while (v.y() < -leny/2.0)
      v[1] += leny;
  }
  if (periodic[2]) {
    while (v.z() > lenz/2.0)
      v[2] -= lenz;
    while (v.z() < -lenz/2.0)
      v[2] += lenz;
  }
  return v;
}

inline double ran(){
  const long unsigned int x = 0;
  static MTRand my_mtrand(x); // always use the same random number generator (for debugging)!
  return my_mtrand.randExc(); // which is the range of [0,1)
}

inline Vector3d ran3(){
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

static void took(const char *name) {
    assert(name); // so it'll count as being used...
    static clock_t last_time = clock();
    clock_t t = clock();
    double peak = peak_memory()/1024.0/1024;
    double seconds = (t-last_time)/double(CLOCKS_PER_SEC);
    if (seconds > 120) {
      printf("%s took %.0f minutes and %g M memory.\n", name, seconds/60, peak);
    } else {
      printf("%s took %g seconds and %g M memory.\n", name, seconds, peak);
    }
    fflush(stdout);
    last_time = t;
}
