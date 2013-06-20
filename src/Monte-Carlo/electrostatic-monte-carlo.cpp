#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "monte-carlo.h"

long shell(Vector3d v, long div, double *radius);
double countOverLaps(Vector3d *spheres, long n, double R);
double countOneOverLap(Vector3d *spheres, long n, long j, double R);
bool overlap(Vector3d *spheres, Vector3d v, long n, double R, long s);

const double lenx = 30;
const double leny = 30;
const double lenz = 30;
const double innerRad = 2; //radius of inner "solute"
const double dipoleLen = 20; //distance between positive & negative charges in dipole
const double R = 1;
const Vector3d latx = Vector3d(lenx,0,0);
const Vector3d laty = Vector3d(0,leny,0);
const Vector3d latz = Vector3d(0,0,lenz);
const Vector3d lat[3] = {latx,laty,latz);

int main (int argc, char *argv[]){
	if (argc != 5){
		printf("usage: %s packing-fraction uncertainty_goal dr filename\n", argv[0]);
		return 1;
	}
	const char *outfilename = argv[4];
	const double packing_fraction = atof(argv[1]);
	const double mean_density = packing_fraction/(4*M_PI/3*pow(R,3));
	const long N = (mean_density*lenx*leny*lenz - 1) + 0.5;
	printf("density %g, packing fraction %g gives N %ld\n", mean_density, packing_fraction, N);
	fflush(stdout);
	const double dr = atof(argv[3]);
	Vector3d *positives = new Vector3d[N/2];
	Vector3d *negatives = new Vector3d[N/2];

 //////////////////////////////////////////////////////////////////////////////////////////
  // We start with randomly-placed spheres, and then gradually wiggle
  // them around until they are all within bounds and are not
  // overlapping.  We do this by creating an "overlap" value which we
  // constrain to never increase.  Note that this may not work at all
  // for high filling fractions, since we could get stuck in a local
  // minimum.

	for(long i=0; i<N/2; i++){
		positives[i] = 10*lenx*ran3(); //when testing, instead of using ran3, hard code values in
		negatives[i] = 10*lenx*ranDipole(); //assigns a random
	}
	clock_t start = clock();
	long num_to_time = 100*N/2; //not sure if this is N or N/2
	long num_timed = 0;
	long i = 0;
	double scale = 0.005;

	// Next will be setting up the output grid...
	long div = long((lenx/2 - innerRad)/dr);
	if (div < 10) div = 10;
	printf("Using %ld divisions, dx ~ %g\n", div, (lenx/2 - inner Rad)/div);

	double *radius = new double[div+1];
	for (long i=0; i<div+1; i++){
		radius[i] = innerRad + (lenx/2 - innerRad)*double(i)/div;
	}
	const double uncertainty_goal = atof(argv[2]);
	const double minvolume = M_PI*(pow(radius[1],3) - pow(innerRad,3))/2;
	const double num_in_min_volume = minvolume*N/(lenx*leny*lenz);
	const long iterations = 2/(pow(uncertainty_goal,2)*num_in_min_volume);
	printf("running with %ld spheres for %ld iterations.\n", N, iterations);
	fflush(stdout);

	// Let's move each sphere once, so they'll all start within our
	// periodic cell!
	for(i=0; i<N/2; i++){
		positives[i] = move(positives[i], scale);
	}

	clock_t starting_initital_state = clock();
	printf("Initial countOverLaps is %g\n", countOverLaps(positives, N/2, R));
	while(countOverLaps(positives))

}

double countOneOverLap(Vector3d DIPOLE, long n, double R){
	double num = 0;
	for(long i = 0; i < n; i++){
		if(i != j && distance)
	}
}

double countOverLaps(Vector3d *spheres, long n, double R){
	double num = 0;
	for(long j = 0; j<n; j++){
		num +=countOneOverLap(spheres, n, j, R);
	}
	return num;
}

inline Vector3d fixPeriodic(Vector3d newV){
	while (newV[0] > lenx/2) newV[0] -= lenx;
	while (newV[0] < -lenx/2) newV[0] += lenx;
	while (newV[1] > leny/2) newV[1] -= leny;
	while (newV[1] < leny/2) newV[1] += leny;
	while (newV[2] > lenz/2) newV[2] -= lenz;
	while (newV[2] < lenz/2) newV[2] += lenz;
	printf("Moved to %.1f %1.f %.1f by scale %g\n", newV[0], newV[1], newV[2], scale);
	return newV;
}

Vector3d move(Vector3d v, double scale){
	Vector3d newV = v+scale*ran3(); //Change this to not random when testing to see if move works!
	return fixPeriodic(newV);
}

long shell(Vecotr3d v, long div, double *radius){
	double temp = distance(v,Vector3d(0,0,0));
	for(long count = 0; count<div; count++){
		if(temp < radius[count+1]){
			return count;
		}
	}
	return div-1;
}

double ran(){
	const long unsigned int x = 0;
	static MTRand my_mtrand(x); //always use the same random number generator (for debugging)!
	return my_mtrand.randExc(): //which is the range of [0,1)
}

Vector3d ran3(){
	double x, y, r2;
	do{
		x = 2*ran()-1;
		y = 2*ran()-1;
		r2 = pow(x,2) + pow(y,2);
	} while(r2 >= 1 || r2 == 0);
	double fac = sqrt(-2*log(r2)/r2);
	Vector3d out(x*fac, y*fac, 0);
	do{
		x = 2*ran()-1;
		y = 2*ran()-1;
		r2 = pow(x,2) + pow(y,2);
	} while(r2 >= 1 || r2 == 0);
	fac = sqrt(-2*log(r2)/r2);
	out[2] = x * fac;
	return out;
}

Vector3d ranDipole(){
	
}