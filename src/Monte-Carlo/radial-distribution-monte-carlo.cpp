#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"
#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
using std::vector;

long shell(Vector3d v, long div, double *radius);
double countOverLaps(Vector3d *spheres, long n, double R);
double countOneOverLap(Vector3d *spheres, long n, long j, double R);
bool overlap(Vector3d *spheres, Vector3d v, long n, double R, long s);

const double lenx = 30;
const double leny = 30;
const double lenz = 30;
const double innerRad = 2;  //of inner spherical "solute"
const double R = 1;
const Vector3d latx = Vector3d(lenx,0,0);
const Vector3d laty = Vector3d(0,leny,0);
const Vector3d latz = Vector3d(0,0,lenz);
const Vector3d lat[3] = {latx,laty,latz};

int main(int argc, char *argv[]){
  if (argc != 5) {
    printf("usage:  %s packing-fraction uncertainty_goal dr filename\n", argv[0]);
    return 1;
  }
  const char *outfilename = argv[4];
  const double packing_fraction = atof(argv[1]);
  const double mean_density = packing_fraction/(4*M_PI/3*R*R*R);
  const long N = (mean_density*lenx*leny*lenz - 1) + 0.5;
  printf("density %g, packing fraction %g gives N %ld\n", mean_density, packing_fraction, N);
  fflush(stdout);
  const double dr = atof(argv[3]);
  Vector3d *spheres = new Vector3d[N];

  //////////////////////////////////////////////////////////////////////////////////////////
  // We start with randomly-placed spheres, and then gradually wiggle
  // them around until they are all within bounds and are not
  // overlapping.  We do this by creating an "overlap" value which we
  // constrain to never increase.  Note that this may not work at all
  // for high filling fractions, since we could get stuck in a local
  // minimum.
  for(long i=0; i<N; i++) {
    spheres[i]=10*lenx*ran3();
  }
  clock_t start = clock();
  long num_to_time = 100*N;
  long num_timed = 0;
  long i = 0;
  double scale = .005;

  // At this stage, we'll set up our output grid...
  long div = long((lenx/2 - innerRad)/dr);
  if (div < 10) div = 10;
  printf("Using %ld divisions, dx ~ %g\n", div, lenx/2/div);

  double *radius = new double[div+1];
  for (long i=0;i<div+1;i++) radius[i] = innerRad + (lenx/2 - innerRad)*double(i)/div;

  const double uncertainty_goal = atof(argv[2]);
  const double minvolume = M_PI*(radius[1]*radius[1]*radius[1] - innerRad*innerRad*innerRad)/2;
  const double num_in_min_volume = minvolume*N/lenx/leny/lenz;
  const long iterations = 2.0/uncertainty_goal/uncertainty_goal/num_in_min_volume;
  printf("running with %ld spheres for %ld iterations.\n", N, iterations);
  fflush(stdout);


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

  // Here we use a hokey heuristic to decide on an average move
  // distance, which is proportional to the mean distance between
  // spheres.
  const double mean_spacing = pow(lenx*leny*lenz/N, 1.0/3);
  if (mean_spacing > 2*R) {
    scale = 2*(mean_spacing - 2*R);
  } else {
    scale = 0.1;
  }
  printf("Using scale of %g\n", scale);
  long count = 0;
  long *shells = new long[div];
  for (long l=0; l<div; l++) shells[l] = 0;

  double *density = new double[div];

  /////////////////////////////////////////////////////////////////////////////
  int hours_now = 1;
  num_to_time = 1000;
  start = clock();
  num_timed = 0;
  double secs_per_iteration = 0;
  long workingmoves=0;

  clock_t output_period = CLOCKS_PER_SEC*60; // start at outputting every minute
  clock_t max_output_period = CLOCKS_PER_SEC*60*60; // top out at one hour interval
  clock_t last_output = clock(); // when we last output data
  for (long j=0; j<iterations; j++){
	  num_timed = num_timed + 1;
    if (num_timed > num_to_time || j == iterations - 1) {
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
      if (now > last_output + output_period || j == iterations - 1) {
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
            printf("Saved data after one hour %ld minutes\n", mins_done);
          } else {
            printf("Saved data after %ld hours, %ld minutes\n", hours_done, mins_done);
          }
          fflush(stdout);
        }
        for(long i=0; i<div; i++){
          const double rmax = radius[i+1];
          const double rmin = radius[i];
          const double num_counted = (j+1)/double(N);
          const double dV = 4/3.*M_PI*(rmax*rmax*rmax - rmin*rmin*rmin);
          density[i]=shells[i]/dV/num_counted*(4*M_PI/3);
        }

        FILE *out = fopen((const char *)outfilename,"w");
        if (out == NULL) {
          printf("Error creating file %s\n", outfilename);
          return 1;
        }
        if (fabs(density[1] - density[0])/density[0] > 0.1/sqrt(shells[0])) {
          // We have enough precision to consider extrapolating to the
          // contact point.  We could do better than this, but this
          // should be good enough when we have solid statistics.
          fprintf(out, "%g\t%g\n", radius[0], 1.5*density[0] - 0.5*density[1]);
        } else {
          // We'll just duplicate the lowest point we have, as if it
          // were the value at contact.
          fprintf(out, "%g\t%g\n", radius[0], density[0]);
        }
        fprintf(out, "%g\t%g\n", 0.5*(radius[0]+radius[1]), density[0]);
        long divtoprint = div;
        divtoprint = div - 1;
        for(long i=1; i<divtoprint; i++){
          fprintf(out, "%g\t%g\n", 0.5*(radius[i]+radius[i+1]), density[i]);
        }
        fflush(stdout);
        fclose(out);
      }
      ///////////////////////////////////////////end of print.dat
    }
	
    // only write out the sphere positions after they've all had a
    // chance to move
    if (workingmoves%N == 0) {
      for (long i=0;i<N;i++) {
        //printf("Sphere at %.1f %.1f %.1f\n", spheres[i][0], spheres[i][1], spheres[i][2]);
        shells[shell(spheres[i], div, radius)]++;
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
      continue;
    }
    spheres[j%N] = temp;
    workingmoves++;
    if (scale < 5 && false) {
      scale = scale*1.02;
      //printf("Increasing scale to %g\n", scale);
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////

  printf("Total number of attempted moves = %ld\n",count);
  printf("Total number of successful moves = %ld\n",workingmoves);
  printf("Acceptance rate = %g\n", workingmoves/double(count));
  fflush(stdout);
  delete[] shells;
  delete[] density;
  delete[] spheres;
}

double countOneOverLap(Vector3d *spheres, long n, long j, double R){
  double num = 0;
  for(long i = 0; i < n; i++){
    if(i != j && distance(spheres[i],spheres[j])<2*R){
      num+=2*R-distance(spheres[i],spheres[j]);
    }
  }
  if(distance(spheres[j],Vector3d(0,0,0))<innerRad){
    num += innerRad - distance(spheres[j],Vector3d(0,0,0));
  }
  if (spheres[j][0] > lenx/2 ){
    num += spheres[j][0]-(lenx/2);
  } else if (spheres[j][0] < -lenx/2){
    num -= spheres[j][0] + (lenx/2);
  }
  if (spheres[j][1] > leny/2 ){
    num += spheres[j][1]-(leny/2);
  } else if (spheres[j][1] < -leny/2){
    num -= spheres[j][1] + (leny/2);
  }
  if (spheres[j][2] > lenz/2 ){
    num += spheres[j][2]-(lenz/2);
  } else if (spheres[j][2] < -lenz/2){
    num -= spheres[j][2] + (lenz/2);
  }
  for(long i = 0; i < n; i++){
    if (i != j) {
      for (long k=0; k<3; k++){
        if (distance(spheres[j],spheres[i]+lat[k]) < 2*R){
          num += 2*R - distance(spheres[j],spheres[i]+lat[k]);
        } else if (distance(spheres[j],spheres[i]-lat[k]) < 2*R) {
          num += 2*R - distance(spheres[j],spheres[i]-lat[k]);
        }
        for (long m=k+1; m<3; m++){
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
  if (distance(v,Vector3d(0,0,0)) < innerRad) return true;
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
    for(long i = 0; i < n; i++) {
      if (i!=s){
        if (distance(v,spheres[i]+lat[k]) < 2*R) return true;
        if (distance(v,spheres[i]-lat[k]) < 2*R) return true;
      }
    }
    for (long m=k+1; m<3; m++){
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
  if (amonborder[0] && amonborder[1] && amonborder[2]){
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
  while (newv[0] > lenx/2) newv[0] -= lenx;
  while (newv[0] < -lenx/2) newv[0] += lenx;
  while (newv[1] > leny/2) newv[1] -= leny;
  while (newv[1] < -leny/2) newv[1] += leny;
  while (newv[2] > lenz/2) newv[2] -= lenz;
  while (newv[2] < -lenz/2) newv[2] += lenz;
  //printf("Moved to %.1f %.1f %.1f by scale %g\n", newv[0], newv[1], newv[2], scale);
  return newv;
}

Vector3d move(Vector3d v,double scale){
  Vector3d newv = v+scale*ran3();
  return fixPeriodic(newv);
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
    if (distance(v,w+lat[k]) < 2*oShell) return true;
    if (distance(v,w-lat[k]) < 2*oShell) return true;
    for (long m=k+1; m<3; m++){
      if (distance(v,w+lat[k]+lat[m]) < 2*oShell) return true;
      if (distance(v,w-lat[k]-lat[m]) < 2*oShell) return true;
      if (distance(v,w+lat[k]-lat[m]) < 2*oShell) return true;
      if (distance(v,w-lat[k]+lat[m]) < 2*oShell) return true;
    }
  }
  if (distance(v,w+latx+laty+latz) < 2*oShell) return true;
  if (distance(v,w+latx+laty-latz) < 2*oShell) return true;
  if (distance(v,w+latx-laty+latz) < 2*oShell) return true;
  if (distance(v,w-latx+laty+latz) < 2*oShell) return true;
  if (distance(v,w-latx-laty+latz) < 2*oShell) return true;
  if (distance(v,w-latx+laty-latz) < 2*oShell) return true;
  if (distance(v,w+latx-laty-latz) < 2*oShell) return true;
  if (distance(v,w-latx-laty-latz) < 2*oShell) return true;
  return false;
}



long shell(Vector3d v, long div, double *radius){
  double temp = distance(v,Vector3d(0,0,0));
  for(long count = 0; count<div; count++){
    if(temp<radius[count+1]) return count;
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
