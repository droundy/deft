#include <stdio.h>
#include <time.h>
#include <Eigen/Core>
#include "utilities.h"
#include "MersenneTwister.h"
USING_PART_OF_NAMESPACE_EIGEN

#include <cassert>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double ran();
Vector3d ran3();

inline double distance(Vector3d v1, Vector3d v2){
  return (v1 - v2).norm();
}

Vector3d move(Vector3d v, double R);
long shell(Vector3d v, long div, double *radius, double *sections);
bool overlap(Vector3d *spheres, Vector3d v, long n, double R, long s);
double calcPressure(Vector3d *spheres, long N, double volume);
double potentialEnergy(Vector3d *spheres, long n, double R);
inline Vector3d fixPeriodic(Vector3d newv);

double kT;
double eps = 1;
double rho = 1;
bool testp = false;
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
const double R = 1;
double R_wall = R;
double length_of_cavity=0;
double R_T = R + R_wall;
double testp_sigma, testp_eps, testp_r;
Vector3d latx = Vector3d(lenx,0,0);
Vector3d laty = Vector3d(0,leny,0);
Vector3d latz = Vector3d(0,0,lenz);
Vector3d lat[3] = {latx,laty,latz};
bool FCC=false;  //if true then face centered cubic latice
bool flat_div = false; //the divisions will be equal and will divide from z wall to z wall
bool WCA = true; //Defaults to a test particle with Lennard-Jones repulsion
bool soft_wall = false;
const double sigma = R*pow(2,5.0/6.0);
double sig_wall = R_T*pow(2,-1.0/6.0);

bool periodic[3] = {false, false, false};
inline double max(double a, double b) { return (a>b)? a : b; }

int main(int argc, char *argv[]){
  if (argc < 6) {
    printf("usage:  %s Nspheres dx uncertainty_goal filename \n there will be more!\n", argv[0]);
    return 1;
  }


  printf("Command line: ");
  for (int i=0;i<argc;i++) {
    printf("%s ", argv[i]);
  }
  printf("\n");
 
  WCA = true;
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
      flat_div = true; // we want flat divisions in the z directoin,
                       // if we have a z wall.  In every other case we
                       // use spherical divisions.
      lenz = atof(argv[a+1]);
      periodic[2] = false;
      maxrad = max(maxrad, lenz);
    } else if (strcmp(argv[a],"softwallz") == 0) {
      has_z_wall = true;
      flat_div = true; // we want flat divisions in the z directoin,
                       // if we have a z wall.  In every other case we
                       // use spherical divisions.
      R_wall = atof(argv[a+1]);
      lenz = 30;
      periodic[2] = false;
      maxrad = max(maxrad, lenz);
      soft_wall = true;
      R_T = R + R_wall;
      sig_wall = R_T*pow(2,-1.0/6.0);
    } else if (strcmp(argv[a],"kT") == 0) {
      kT = atof(argv[a+1]);
      printf("kT is %g\n",kT);
    } else if (strcmp(argv[a],"TestP")==0){
        testp_sigma = atof(argv[a+1]);
        testp = true;
        testp_r = pow(2,-6/2)*testp_sigma;
    } else if (strcmp(argv[a],"testp_eps")==0){
	testp_eps = atof(argv[a+1]);
        printf("Testp is %d, sigma is %f, eps is %f, last thing is %s\n",testp, testp_sigma,testp_eps,argv[a+1]);
    } else if (strcmp(argv[a],"potential") == 0){
      if (strcmp(argv[a+1],"wca") == 0) {
        WCA = true;
        printf("Using Weeks-Chandler-Andersen potential\n");
      } else {
        printf("Using Weeks-Chandler-Andersen potential by default\n");
      }
    } else if (strcmp(argv[a],"fcc")==0){
	FCC=true;
	length_of_cavity = atof(argv[a+1]);
        lenx = leny = lenz = length_of_cavity;
        rad = maxrad = lenx/2;
        periodic[0] = true;
        periodic[1] = true;
        periodic[2] = true;
	printf("FCC is true length of cavity is %g\n", length_of_cavity);
    }  else {
      printf("Bad argumendat:  %s\n%d\n", argv[a],(strcmp(argv[a],"fcc")));
	printf("%d" ,(strcmp(argv[a],"fcc")));
	return 1;
      }
  }
  fflush(stdout);

  //printf("flatdiv = %s\n", flat_div ? "true" : "false");
  //printf("outerSphere = %s\n", spherical_outer_wall ? "true" : "false");
  //printf("innerSphere = %s\n", spherical_inner_wall ? "true" : "false");
  latx = Vector3d(lenx,0,0);
  laty = Vector3d(0,leny,0);
  latz = Vector3d(0,0,lenz);
  lat[0] = latx;
  lat[1] = laty;
  lat[2] = latz;
  double volume = lenx*leny*lenz;
  if (spherical_inner_wall){
    volume = lenx*leny*lenz-(4*M_PI*(1/3))*(innerRad*innerRad*innerRad);
    //Possible error in these volume calculations.  With
    //spherical_inner_wall: volume should be
    //lenx*leny*lenz-(4*M_PI*(1/3))*(innerRad*innerRad*innerRad) -Sam
  } else if (spherical_outer_wall) {
    volume = 4*M_PI*(1/3)*(rad*rad*rad);
    //Possible error in these volume calculations. With
    //spherical_outer_wall, volume should be
    //4*M_PI*(1/3)*(rad*rad*rad)... Also add need to add a volume for when both
    //are true.  -Sam
  } 

  const char *outfilename = argv[4];

  long N = atol(argv[1]);
  if (FCC == true){
    N = 4*7*7*7;
  }
  printf("firs ttime N is %ld\n",N);
  const double dx_goal = atof(argv[2]);
  const double uncertainty_goal = atof(argv[3]);
  long density_saved_count = 0;

  Vector3d *spheres = new Vector3d[N];
  const int num_frames = 100;
  int which_frame = 0;
  Vector3d *movieData = new Vector3d[num_frames*N];

  if (uncertainty_goal < 1e-12 || uncertainty_goal > 1.0) {
    printf("Crazy uncertainty goal:  %s\n", argv[1]);
    return 1;
  }
  printf("running with %ld spheres.\n", N);

  //////////////////////////////////////////////////////////////////////////////////////////
  // We start with randomly-placed spheres, and then gradually wiggle
  // them around until they are all within bounds and are not
  // overlapping.  We do this by creating an "overlap" value which we
  // constrain to never increase.  Note that this may not work at all
  // for high filling fractions, since we could get stuck in a local
  // minimum.
  // start of initialization
  clock_t start = clock();
  long num_to_time = 100*N;
  long num_timed = 0;
  const long iterations_per_pressure_check = 10*N;
  double scale = .1;
  // Here we use a hokey heuristic to decide on an average move
  // distance, which is proportional to the mean distance between
  // spheres.
  const double mean_spacing = pow(rad*rad*rad/N, 1.0/3);
  if (mean_spacing > 2*R) {
    scale = 2*(mean_spacing - 2*R);
  }
  printf("Using scale of %g (mean spacing = %g)\n", scale, mean_spacing);
  long iterations = 0;
  
  // Let's move each sphere once, so they'll all start within our
  // periodic cell!

  if (FCC == true) {
    printf("testkljdfjak\n");
    iterations = 9999999999;
    //    int cubeLength = pow(double(N/4.0),1.0/3.0) + .5;
    //    int numberOfCubes = ceil(N/4);
    int x = 0;
    int y = 0;
    int z = 0;
    int x1 = 0;
    int y1 = 0;
    int z1 = 0;
    int round = 0;
    for(long k=0;k<N;k++){
      spheres[k][0]=x*(length_of_cavity/7.0) + x1*(1.0/2.0)*(length_of_cavity)*(1.0/7.0);
      spheres[k][1]=y*(length_of_cavity/7.0) + y1*(1.0/2.0)*(length_of_cavity)*(1.0/7.0);
      spheres[k][2]=z*(length_of_cavity/7.0) + z1*(1.0/2.0)*(length_of_cavity)*(1.0/7.0);
      spheres[k] = fixPeriodic(spheres[k]); // move sphere into the periodic cell
      x=x+1;
      if (x >= 7){//pow(double(numberOfCubes),(1.0/3.0))){
	x=0;
	y=y+1;
      }
      if (y >= 7){//pow(double(numberOfCubes),(1.0/3.0))){
	y=0;
	z=z+1;
      }
      if (z >= 7){
	round = round + 1;
	x=0;
	y=0;
	z=0;
	if (round <= 1){
	  x1=1;
	  y1=1;
	  z1=0;
	  
	}else if (round == 2){
	  x1=1;
	  y1=0;
	  z1=1;
	}else if (round == 3){
	  x1=0;
	  y1=1;
	  z1=1;	  
	}
	//	assert(round<=3 || k == N-1);
      }
    }

    char *gfilename = new char[1000];
    sprintf(gfilename, "%s.positions", outfilename);
    FILE* out = fopen(gfilename, "w");
    if (out == NULL) {
      printf("Error creating file %s\n", outfilename);
      return 1;
    }
    for(long i=0;i<N;i++){
      fprintf(out, "%g\t%g\t%g\n", spheres[i][0], spheres[i][1], spheres[i][2]);
    }	
    fclose(out);
    fflush(stdout);
  
  }else {
    
    for(long i=0; i<N; i++) { 
      spheres[i]=rad*ran3();
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

    for (int i=0;i<N;i++) spheres[i] = move(spheres[i], scale);
    {
      // First we'll run the simulation a while to get to a decent
      // starting point...
      long paranoia = 50; // How careful must we be to wait a long time while initializing.
      double oldPE = potentialEnergy(spheres,N,R);
      double newPE = oldPE, olderPE;
      long counter = 0;
      long iters_initializing = 0, successes_initializing = 0;
      do {
	counter = counter + 1;
	olderPE = oldPE;
	oldPE = newPE;
	for (int j=0;j<paranoia*iterations_per_pressure_check/N;j++) {
	  for (int i=0;i<N;i++) {
	    Vector3d temp = move(spheres[i],scale);
	    if(!overlap(spheres, temp, N, R, i)){
	      spheres[i]=temp;
	      successes_initializing++;
	    }
	    iters_initializing++;
	  }
	}
	newPE = potentialEnergy(spheres,N,R);
	const double pressure = calcPressure(spheres, N, volume);
	if (counter%1 == 0) {
	  printf("Potential energy is %g and pressure is %g (my paranoia = %ld)\n",
		 newPE, pressure, paranoia);
	  fflush(stdout);
	}
	paranoia = long(paranoia*1.3); // get a little more paranoid each time...
	
	// We keep going until the energy increases twice.  This is a
	// heuristic to try to avoid starting taking data before the
	// system has reached equilibrium.
      } while (olderPE > oldPE || oldPE > newPE);
      printf("olderPE %g  oldPE %g  newPE %g\n", olderPE, oldPE, newPE);
      
      // At this stage, iters_initializing is the number of iterations
      // we spent trying to get the system into a decent starting state.
      // It seems risky to not run the simulation much more than we
      // needed just to get it to a reasonable starting state.  So here
      // we set the number of iterations accordingly.
      iterations = 4*(iters_initializing + iterations_per_pressure_check);
      printf("running with an extra %ld iterations.\n", iterations);
      printf("success rate initializing was %.0f\n", 100*double(successes_initializing)/iters_initializing);
    }
    if (dx_goal == 0) {
      // This indicates that we are not interested in density variation
      // at all!
      iterations += 1.0/(N*uncertainty_goal*uncertainty_goal);
    } else {
      // The following should give us approaximately uncertainty_goal
      // uncertainty in the density in each bin, assuming all the bins
      // have an equal volume.
      iterations += 1.0/((dx_goal/maxrad)*uncertainty_goal*uncertainty_goal);
    }
    printf("Running for a total of %ld iterations.\n", iterations);
  } //end of initialization
  long div = maxrad/dx_goal;
  if (div < 10) div = 10;
  printf("Using %ld divisions, dx ~ %g\n", div, maxrad/div);
  printf("test abc\n");
  fflush(stdout);
  double *radius = new double[div+1];
  double *sections = new double [div+1];
  double *shellsRadius = new double[div+1];
  double FCC_div = pow(div,.5);
  FCC_div = ceil(FCC_div);
  long FCC_div2 = pow(FCC_div,2)+2*FCC_div;
  double **FCC_shells;
  FCC_shells = new double*[FCC_div2+1];
  printf("fcc_div is %g fcc_div2 is %ld\n",FCC_div,FCC_div2);
  for (long k=0;k<FCC_div2+1;k++){
    FCC_shells[k] = new double[3];
  }
 //divisions will be a fcc_div x fcc_div 2d map
  printf("test def\n");
  long note = 0;
  long note1 = 0;
  if (FCC){
    for (long s=0;s<FCC_div2+1;s++){
      FCC_shells[s][2]=0;
      FCC_shells[s][0] = note1*length_of_cavity*(1/FCC_div);
      FCC_shells[s][1] = note*length_of_cavity*(1/FCC_div);
      note1=note1+1;
      if (s>(note+1)*(FCC_div+1)){
	//printf("true s is %ld note + 1 is %ld note1 is %ld\n", s, note+1, note1); 
	note = note+1;
	note1 = 0;
      }
      
    }
  }
  if (flat_div){
    double size = lenz/div;
    for (long s=0; s<div+1; s++){
      sections[s] = size*s - lenz/2.0;
    }
  }
  if (spherical_inner_wall || testp){
    double size = rad/div;
    for (long s=0; s<div+1; s++){
      radius[s] = size*s;
      shellsRadius[s]=radius[s];
    }
  } else {
    double size = rad/div;
    const double w = 1.0/(1 + dx_goal*div);
    for (long l=0;l<div+1;l++) {
      // make each bin have about the same volume
      shellsRadius[l]=size*l;
      //printf("shells radius is %g\n", shellsRadius[l]);
      radius[l] = w*rad*pow(double(l)/(div), 1.0/3.0) + (1-w)*rad*double(l)/div;
    }
  }
  
  double shellTotalVolume = M_PI*(4/3)*(shellsRadius[div]*shellsRadius[div]*shellsRadius[div]);

  printf("shellsRadius[dv+1] is is is  %g\n", shellsRadius[div]);
  printf("shell total volume %g\n", shellTotalVolume);

  long count = 0;
  long *shells = new long[div];
  long *radial_distributon_histogram = new long[div];

  double *density = new double[div];
  double *radial_distribution = new double[div+1];
  for (long l=0; l<div; l++) {
    shells[l] = 0;
    radial_distributon_histogram[l]=0;
  }

  num_to_time = 5000;
  start = clock();
  num_timed = 0;
  double secs_per_iteration = 0;
  long workingmoves=0;
  long num_pressures_in_sum=0;
  double pressure_sum = 0;

	
  clock_t output_period = CLOCKS_PER_SEC*60; // start at outputting every minute
  clock_t max_output_period = clock_t(CLOCKS_PER_SEC)*60*30; // top out at one hour interval
  clock_t last_output = clock(); // when we last output data
  for (int i=0;i<N;i++) {
    movieData[i] = spheres[i];
    // printf("lenx: %g leny: %g\t\tx: %g y: %g\t\tx': %g y': %g\n",
    //        lenx, leny, movieData[i][0], movieData[i][1],
    //        fixPeriodic(movieData[i])[0], fixPeriodic(movieData[i])[1]);
    assert(fabs(movieData[i][0]) <= lenx/2);
    assert(fabs(movieData[i][1]) <= leny/2);
  }
  which_frame = N;
  for (long j=0; j<iterations; j++){
    num_timed = num_timed + 1;
    if (num_timed > num_to_time || j==(iterations-1)) {
      num_timed = 0;
      fflush(stdout);
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
      if (now - last_output > output_period || j == iterations-1) {
        last_output = now;
        const double increase_ratio = 1.2;
        if (output_period < max_output_period/increase_ratio) {
          output_period *= increase_ratio;
        } else {
          output_period = max_output_period;
        }
        {
          double secs_done = double(now)/CLOCKS_PER_SEC;
          long mins_done = secs_done / 60;
          long hours_done = mins_done / 60;
          mins_done = mins_done % 60;
          long percent_complete = j/iterations*100;
          if (hours_done > 50) {
            printf("Saved data after %ld hours\n", hours_done);
          } else if (mins_done < 1) {
            printf("Saved data after %.1f seconds\n", secs_done);
          } else if (hours_done < 1) {
            printf("Saved data after %ld minutes\n", mins_done);
          } else if (hours_done < 2) {
            printf("Saved data after %ld minutes\n", mins_done);
          } else {
            printf("Saved data after %ld hours, %ld minutes\n", hours_done, mins_done);
          }
          printf("Current/Max Iterations: %ld/%ld\n", j, iterations);
          printf("Percent Complete: %ld  (acceptance rate = %g)\n", percent_complete, workingmoves/double(count));
        }
        for (long i=0; i<div; i++) {
          const double Vi = (shellsRadius[i+1]*shellsRadius[i+1]*shellsRadius[i+1]-shellsRadius[i]*shellsRadius[i]*shellsRadius[i])*(4/3.*M_PI);
          radial_distribution[i] = radial_distributon_histogram[i]*volume/(Vi*density_saved_count*N*(N-1));
        }
        if (!flat_div){
          for(long i=0; i<div; i++){
            double rmax = radius[i+1];
            double rmin = radius[i];
            const double dvolume = 4/3.*M_PI*rmax*rmax*rmax - 4/3.*M_PI*rmin*rmin*rmin;
            density[i]=shells[i]/dvolume/density_saved_count;
          }
        } else {
          const double dvolume = lenx*leny*lenz/div;
          for(long i=0; i<div; i++){
            density[i]=shells[i]/dvolume/density_saved_count;
          }
        }

        FILE *out = fopen(outfilename,"w");
        if (out == NULL) {
          printf("Error creating file %s\n", outfilename);
          return 1;
        }
        if (flat_div){
          fprintf(out, "%g\t%g\n", 0.5*(sections[0]+sections[1]), density[0]);
        } else if (spherical_inner_wall || testp) {
          fprintf(out, "%g\t%g\n", radius[0], 0.0);
          fprintf(out, "%g\t%g\n", 0.5*(radius[0]+radius[1]), density[0]);
        } else {
          fprintf(out, "%g\t%g\n" , 0.0, density[0]);
        }
        long divtoprint = div;
        if (!spherical_outer_wall) divtoprint = div - 1;
        if (!flat_div) {
          for(long i=1; i<divtoprint; i++){
            fprintf(out, "%g\t%g\n", 0.5*(radius[i]+radius[i+1]), density[i]);
          }
        } else {
          for(long i=1; i<div; i++){
            fprintf(out, "%g\t%g\n", 0.5*(sections[i]+sections[i+1]), density[i]);
          }
        }
        fclose(out);



        if (periodic[0] && periodic[1] && periodic[2]) {
          char *gfilename = new char[1000];
          sprintf(gfilename, "%s.gradial", outfilename);
          out = fopen(gfilename, "w");
          if (out == NULL) {
            printf("Error creating file %s\n", gfilename);
            return 1;
          }
          delete[] gfilename;
          fprintf(out, "%g\t%g\n", 0.0, radial_distribution[0]);
          for (long i=1; i<div; i++) {
            fprintf(out, "%g\t%g\n", 0.5*(shellsRadius[i-1]+shellsRadius[i]), radial_distribution[i]);
          }
          fclose(out);
        }
	
	if (FCC) {
	  
          char *gfilename = new char[1000];
          sprintf(gfilename, "%s.fcc", outfilename);
          out = fopen(gfilename, "w");
          if (out == NULL) {
            printf("Error creating file %s\n", gfilename);
            return 1;
          }
          delete[] gfilename;
	  //          fprintf(out, "%g\t%g\n", 0.0, radial_distribution[0]);
	  for (long i=0; i<FCC_div2+1; i++) {
	    fprintf(out, "%g\t%g\t%g\n",FCC_shells[i][0], FCC_shells[i][1], FCC_shells[i][2]);
	  }
          fclose(out);
	  printf("SAVED FCC\n");
        }

        fflush(stdout);
        if (num_pressures_in_sum > 0) {
          char *pressureFileName = new char[10000];
          sprintf(pressureFileName, "%s.prs", outfilename);
          FILE *pressureFile = fopen(pressureFileName, "w");
          if (pressureFile == NULL) {
            printf("Error creating file %s\n", pressureFileName);
            return 1;
          }
          delete[] pressureFileName;
          fprintf(pressureFile, "%g\n", pressure_sum/num_pressures_in_sum);
          fclose(pressureFile);
        } else {
          printf("No pressure data yet after %ld iterations, compared with %ld\n",
                 j, iterations_per_pressure_check);
          fflush(stdout);
        }
        
        fflush(stdout);
      }
      ///////////////////////////////////////////end of print.dat
    }
    
    Vector3d temp = move(spheres[j%N],scale);
    count++;
    if(!overlap(spheres, temp, N, R, j%N)){
      spheres[j%N] = temp;
      workingmoves++;
    }
    if (which_frame < N*num_frames) {
      assert(which_frame%N == j%N);
      movieData[which_frame] = spheres[which_frame%N];
      which_frame++;
    } else if (movieData) {
      const double frame_duration = 500.0; // milliseconds
      printf("\n\nworking on movie!!!\n\n\n");
      const double scale = 10;
      char *movie_fname = (char *)malloc(4000);
      sprintf(movie_fname, "%s.svg", outfilename);
      FILE *movie = fopen(movie_fname, "w");
      free(movie_fname);
      fprintf(movie, "<svg version=\"1.1\" width=\"%g\" height=\"%g\" xmlns=\"http://www.w3.org/2000/svg\">\n",
              lenx*scale, leny*scale);
      for (int s=0;s<N;s++) {
        fprintf(movie,
                "  <circle cx=\"%g\" cy=\"%g\" r=\"%g\" fill=\"none\" stroke=\"black\" stroke-width=\"%g\">\n",
                (lenx/2+movieData[s][0])*scale, (leny/2+movieData[s][1])*scale,
                R*scale, 0.05*R*scale);

        fprintf(movie,
                "    <animate attributeName=\"cx\" dur=\"%gms\" repeatCount=\"indefinite\" values=\"%g",
                frame_duration*num_frames, (lenx/2+movieData[s][0])*scale);
        for (int f=1;f<num_frames;f++) {
          fprintf(movie, "; %g", (lenx/2+movieData[s+f*N][0])*scale);
        }
        fprintf(movie, "; %g\" keyTimes=\"0",
                (lenx/2+movieData[s+(num_frames-1)*N][0])*scale);
        for (int f=1;f<num_frames;f++) {
          fprintf(movie, "; %g", (1+s+(f-1)*N)/double(num_frames*N));
        }
        fprintf(movie, "; 1\"/>\n");

        fprintf(movie,
                "    <animate attributeName=\"cy\" dur=\"%gms\" repeatCount=\"indefinite\" values=\"%g",
                frame_duration*num_frames, (leny/2+movieData[s][1])*scale);
        for (int f=1;f<num_frames;f++) {
          fprintf(movie, "; %g", (leny/2+movieData[s+f*N][1])*scale);
        }
        fprintf(movie, "; %g\" keyTimes=\"0",
                (lenx/2+movieData[s+(num_frames-1)*N][1])*scale);
        for (int f=1;f<num_frames;f++) {
          fprintf(movie, "; %g", (1+s+(f-1)*N)/double(num_frames*N));
        }
        fprintf(movie, "; 1\"/>\n");

        fprintf(movie,
                "    <animate attributeName=\"opacity\" dur=\"%gms\" repeatCount=\"indefinite\" values=\"",
                frame_duration*num_frames);
        bool last_frame_moved = true;
        bool penultimate_frame_moved = false;
        for (int f=1;f<num_frames-1;f++) {
          const double dx = movieData[s+f*N][0]-movieData[s+(f+1)*N][0];
          const double dy = movieData[s+f*N][1]-movieData[s+(f+1)*N][1];
          if (last_frame_moved || penultimate_frame_moved) {
            fprintf(movie, "0; ");
            penultimate_frame_moved = last_frame_moved;
            last_frame_moved = false;
          } else if (fabs(dx) > lenx/2 || fabs(dy) > leny/2) {
            penultimate_frame_moved = last_frame_moved;
            fprintf(movie, "0; ");
          } else {
            penultimate_frame_moved = last_frame_moved;
            fprintf(movie, "1; ");
          }
          if (fabs(dx) > lenx/2 || fabs(dy) > leny/2) {
            last_frame_moved = true;
          }
        }
        fprintf(movie, " 0; 0\"/>\n");

        fprintf(movie, "  </circle>\n");
      }
      fprintf(movie, "</svg>\n");
      fclose(movie);
      delete[] movieData;
      movieData = 0;
    }

    // only write out the sphere positions after they've all had a
    // chance to move
    const long first_data_iterations = 20*iterations_per_pressure_check;
    if (j%N == 0 && workingmoves > first_data_iterations) {
      density_saved_count++;
      for (long s=0;s<N;s++) {
        shells[shell(spheres[s], div, radius, sections)]++;
	if (FCC){
	  for(long n=0;n<FCC_div2;n++){
	    if (FCC_shells[n][0] <= spheres[s][0] &&
                FCC_shells[n+1][0]  >= spheres[s][0] &&
                FCC_shells[n][1] <= spheres[s][1] &&
                FCC_shells[n+1][1]  >= spheres[s][1] ){
	      FCC_shells[n][2] = FCC_shells[n][2]+1;
	      n=FCC_div2+2;
            }
          }
        }
	for (long i=0; i<N; i++){
	  for (long k=0; k<div; k++) {
	    Vector3d vri = spheres[i]-spheres[s];
	    vri = fixPeriodic(vri);
	    const double ri = distance(vri,Vector3d(0,0,0));
	    if (ri < shellsRadius[k+1] && ri > shellsRadius[k] && s != i) {
	      radial_distributon_histogram[k]++;
	    }
	  }
	}
      }
    }
    //  printf("Iterations per pressure %ld working moves is %ld\n",iterations_per_pressure_check, workingmoves);
    if (j%iterations_per_pressure_check == 0 && workingmoves > first_data_iterations) {
      double newpress = calcPressure(spheres, N, volume);
      //const double excpress = newpress - (N/volume)*kT; // difference from ideal gas pressure
      if (newpress != newpress) {
        printf("Got NaN trouble.\n");
        exit(1);
      }
      //printf("Pressure is %g (excess pressure: %g), energy is %g\n",
      //       newpress, excpress, potentialEnergy(spheres,N,R));
      //  for(long i=0; i<div; i++){
        //printf("Number of spheres in division %ld = %ld\n", i+1, shells[i]);
      //}
      pressure_sum += newpress;
      num_pressures_in_sum += 1;
    }
    
    if(j % (first_data_iterations/100)==0 && j < first_data_iterations && j != 0){
      double secs_to_go = secs_per_iteration*(first_data_iterations - j);
      long mins_to_go = secs_to_go / 60;
      long hours_to_go = mins_to_go / 60;
      mins_to_go = mins_to_go % 60;
      if (hours_to_go > 5) {
        printf("%.0f%% to first data... (%ld hours to go)\n",j/(first_data_iterations*1.0)*100, hours_to_go);
      } else if (mins_to_go < 1) {
        printf("%.0f%% to first data... (%.1f seconds to go)\n",j/(first_data_iterations*1.0)*100, secs_to_go);
      } else if (hours_to_go < 1) {
        printf("%.0f%% to first data... (%ld minutes to go)\n",j/(first_data_iterations*1.0)*100, mins_to_go);
      } else if (hours_to_go < 2) {
        printf("%.0f%% to first data... (1 hour, %ld minutes to go)\n",j/(first_data_iterations*1.0)*100, mins_to_go);
      } else {
        printf("%.0f%% to first data... (%ld hours, %ld minutes to go)\n",
               j/(first_data_iterations*1.0)*100, hours_to_go, mins_to_go);
      }
      fflush(stdout);
    }
  }
  ///////////////////////////////////(4*eps*(pow(sigma/r,12) - pow(sigma/r,6)) + eps)/36///////////////////////////////////////////////////////
  for(long i=0; i<div; i++){
    printf("Number of spheres in division %ld = %ld\n  radius is %g\n", i+1, shells[i], radius[i]);
  }

  printf("Total number of attempted moves = %ld\n",count);
  printf("Total number of successful moves = %ld\n",workingmoves);
  printf("Acceptance rate = %g\n", workingmoves/double(count));
  //delete[] shells;
  fflush(stdout);
  //delete[] density;
  fflush(stdout);
  delete[] spheres;

  fflush(stdout);
}

inline double sqr(double x) {
  return x*x;
}

inline double potential(double r) {
  if (r >= 2*R) return 0;
  //eps is defined to give the same curvature at r=2R for both potentials
  return (4*eps*(pow(sigma/r,12) - pow(sigma/r,6)) + eps);
}

inline double soft_wall_potential(double z) {
  if (z < 0) {
    z = lenz/2 + z;
  } else {
    z = lenz/2 - z;
  }
  if (z >  R_T) return 0;
  return 2*M_PI*rho*eps*((pow(z,3) - pow(R_T,3))/6 + 2*pow(sig_wall,12)*(1/pow(z,9)-1/pow(R_T,9))/45 + (R_T-z)*(R_T*R_T/2 + sig_wall*sig_wall*pow(sig_wall/R_T,4) - 2*sig_wall*sig_wall*pow(sig_wall/R_T,10)/5) + pow(sig_wall,6)*(1/pow(R_T,3)-1/pow(z,3))/3);
}

bool overlap(Vector3d *spheres, Vector3d v, long n, double R, long s){
  double energyNew = 0.0;
  double energyOld = 0.0;
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
  bool wasonborder[3] = {
    fabs(spheres[s][0]) + 2*R >= lenx/2,
    fabs(spheres[s][1]) + 2*R >= leny/2,
    fabs(spheres[s][2]) + 2*R >= lenz/2
  };

  // Energy before potential move
  if (testp){
    double r0 = v.norm();	
    energyOld += -(4*testp_eps*(pow(testp_sigma/r0,12) - pow(testp_sigma/r0,6)));
  }
  if (soft_wall) { energyOld += soft_wall_potential(spheres[s][2]); }
  for(long i = 0; i < n; i++){
    if (i != s) energyOld += potential(distance(spheres[i],spheres[s]));
  }
  
  for (long k=0; k<3; k++) {
    if (periodic[k] && wasonborder[k]) {
      for(long i = 0; i < n; i++) {
        if (i != s) {
          energyOld += potential(distance(spheres[s],spheres[i]+lat[k]));
          energyOld += potential(distance(spheres[s],spheres[i]-lat[k]));
        }
      }
      for (long m=k+1; m<3; m++){
        if (periodic[m] && wasonborder[m]){
          for(long i = 0; i < n; i++) {
            if (i != s) {
              energyOld += potential(distance(spheres[s],spheres[i]+lat[k]+lat[m]));
              energyOld += potential(distance(spheres[s],spheres[i]-lat[k]-lat[m]));
              energyOld += potential(distance(spheres[s],spheres[i]+lat[k]-lat[m]));
              energyOld += potential(distance(spheres[s],spheres[i]-lat[k]+lat[m]));
            }
          }
        }
      }
    }
  }
  if (periodic[0] && periodic[1] && periodic[2]
      && wasonborder[0] && wasonborder[1] && wasonborder[2]){
    for(long i = 0; i < n; i++) {
      if (i != s) {
        energyOld += potential(distance(spheres[s],spheres[i]+latx+laty+latz));
        energyOld += potential(distance(spheres[s],spheres[i]+latx+laty-latz));
        energyOld += potential(distance(spheres[s],spheres[i]+latx-laty+latz));
        energyOld += potential(distance(spheres[s],spheres[i]-latx+laty+latz));
        energyOld += potential(distance(spheres[s],spheres[i]-latx-laty+latz));
        energyOld += potential(distance(spheres[s],spheres[i]-latx+laty-latz));
        energyOld += potential(distance(spheres[s],spheres[i]+latx-laty-latz));
        energyOld += potential(distance(spheres[s],spheres[i]-latx-laty-latz));
      }
    }
  }
  // Energy after potential move
  if (testp){  
    double r1 = spheres[s].norm();	
    energyNew += -(4*testp_eps*(pow(testp_sigma/r1,12) - pow(testp_sigma/r1,6)));
  }
  if (soft_wall) { energyNew += soft_wall_potential(v[2]); }
  for(long i = 0; i < n; i++) {
    if (i != s) energyNew += potential(distance(spheres[i],v));
  }

  for (long k=0; k<3; k++) {
    if (periodic[k] && amonborder[k]) {
      for(long i = 0; i < n; i++) {
        if (i != s) {
          energyNew += potential(distance(v,spheres[i]+lat[k]));
          energyNew += potential(distance(v,spheres[i]-lat[k]));
        }
      }
      for (long m=k+1; m<3; m++){
        if (periodic[m] && amonborder[m]){
          for(long i = 0; i < n; i++) {
            if (i != s) {
              energyNew += potential(distance(v,spheres[i]+lat[k]+lat[m]));
              energyNew += potential(distance(v,spheres[i]-lat[k]-lat[m]));
              energyNew += potential(distance(v,spheres[i]+lat[k]-lat[m]));
              energyNew += potential(distance(v,spheres[i]-lat[k]+lat[m]));
            }
          }
        }
      }
    }
  }
  if (periodic[0] && periodic[1] && periodic[2]
      && amonborder[0] && amonborder[1] && amonborder[2]){
    for(long i = 0; i < n; i++) {
      if (i != s) {
        energyNew += potential(distance(v,spheres[i]+latx+laty+latz));
        energyNew += potential(distance(v,spheres[i]+latx+laty-latz));
        energyNew += potential(distance(v,spheres[i]+latx-laty+latz));
        energyNew += potential(distance(v,spheres[i]-latx+laty+latz));
        energyNew += potential(distance(v,spheres[i]-latx-laty+latz));
        energyNew += potential(distance(v,spheres[i]-latx+laty-latz));
        energyNew += potential(distance(v,spheres[i]+latx-laty-latz));
        energyNew += potential(distance(v,spheres[i]-latx-laty-latz));
      }
    }
  }
  double probabilityOfChange = exp((energyNew-energyOld)/-kT);
  double doesItChange = ran();
  if (doesItChange <= probabilityOfChange) return false;
  else return true;
  return false;
}


double potentialEnergy(Vector3d *spheres, long n, double R){
  double potEnergy = 0.0;
  if (testp){
    for (long b=0;b<n;b++){
      double magnitude = spheres[b].norm();	
      potEnergy += -(4*testp_eps*(pow(testp_sigma/magnitude,12) - pow(testp_sigma/magnitude,6)));
    }
  }
  for (long s=0; s<n; s++){
    bool amonborder[3] = {
      fabs(spheres[s][0]) + 2*R >= lenx/2,
      fabs(spheres[s][1]) + 2*R >= leny/2,
      fabs(spheres[s][2]) + 2*R >= lenz/2
    };
    if (soft_wall) { potEnergy += soft_wall_potential(spheres[s][2]); }
    for(long i = s+1; i < n; i++){
      potEnergy += potential(distance(spheres[i],spheres[s]));
    }
    for (long k=0; k<3; k++) {
      if (periodic[k] && amonborder[k]) {
        for(long i = s+1; i < n; i++) {
          potEnergy += potential(distance(spheres[s],spheres[i]+lat[k]));
          potEnergy += potential(distance(spheres[s],spheres[i]-lat[k]));
        }
        for (long m=k+1; m<3; m++){
          if (periodic[m] && amonborder[m]){
            for(long i = s+1; i < n; i++) {
              potEnergy += potential(distance(spheres[s],spheres[i]+lat[k]+lat[m]));
              potEnergy += potential(distance(spheres[s],spheres[i]-lat[k]-lat[m]));
              potEnergy += potential(distance(spheres[s],spheres[i]+lat[k]-lat[m]));
              potEnergy += potential(distance(spheres[s],spheres[i]-lat[k]+lat[m]));
            }
          }
        }
      }
    }
    if (periodic[0] && periodic[1] && periodic[2]
        && amonborder[0] && amonborder[1] && amonborder[2]){
      for(long i = s+1; i < n; i++) {
        potEnergy += potential(distance(spheres[s],spheres[i]+latx+laty+latz));
        potEnergy += potential(distance(spheres[s],spheres[i]+latx+laty-latz));
        potEnergy += potential(distance(spheres[s],spheres[i]+latx-laty+latz));
        potEnergy += potential(distance(spheres[s],spheres[i]-latx+laty+latz));
        potEnergy += potential(distance(spheres[s],spheres[i]-latx-laty+latz));
        potEnergy += potential(distance(spheres[s],spheres[i]-latx+laty-latz));
        potEnergy += potential(distance(spheres[s],spheres[i]+latx-laty-latz));
        potEnergy += potential(distance(spheres[s],spheres[i]-latx-laty-latz));
      }
    }
  }
  return potEnergy;
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

Vector3d move(Vector3d v,double scale){
  Vector3d newv = v+scale*ran3();
  return fixPeriodic(newv);
}

inline double force_times_distance(double rij) {
  if (rij > 2*R) return 0;
  //these forces are negative for repulsive forces. no particular reason
  return -4*eps*(12*pow(sigma/rij,12) - 6*pow(sigma/rij,6));
}

double calcPressure(Vector3d *spheres, long N, double volume){
  double total_force_times_distance = 0;
  for (long s=0; s<N; s++){
    Vector3d v = spheres[s];
    bool amonborder[3] = {
      fabs(v[0]) + 2*R >= lenx/2,
      fabs(v[1]) + 2*R >= leny/2,
      fabs(v[2]) + 2*R >= lenz/2
      };
    
    for(long i = s+1; i < N; i++)
      total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]));

    for (long k=0; k<3; k++) {
      if (periodic[k] && amonborder[k]) {
        for(long i = s+1; i < N; i++) {
            total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]+lat[k]));
            total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]-lat[k]));
        }
        for (long m=k+1; m<3; m++){
          if (periodic[m] && amonborder[m]){
            for(long i = s+1; i < N; i++) {
              total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]+lat[k]+lat[m]));
              total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]-lat[k]-lat[m]));
              total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]+lat[k]-lat[m]));
              total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]-lat[k]+lat[m]));
            }
          }
        }
      }
    }
    if (periodic[0] && periodic[1] && periodic[2] && amonborder[0] && amonborder[1] && amonborder[2]){
      for(long i = s+1; i < N; i++) {
        total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]+latx+laty+latz));
        total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]+latx+laty-latz));
        total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]+latx-laty+latz));
        total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]-latx+laty+latz));
        total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]-latx-laty+latz));
        total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]-latx+laty-latz));
        total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]+latx-laty-latz));
        total_force_times_distance += force_times_distance(distance(spheres[s],spheres[i]-latx-laty-latz));
      }
    }
  }
  //double pressureValue = (N/volume)*kT - (2*M_PI/3)*(1/(6*volume))*total_force_times_distance*total_force_times_distance*(-2*eps/(2*R));
  return (N/volume)*kT - (1/(3*volume))*total_force_times_distance;
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
