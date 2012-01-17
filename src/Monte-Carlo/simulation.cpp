#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"


//CORNERS!!!!!!!!



//Input "1" for large spherical cavity of radius "rad"
//Input "2" for two dimensions boundary conditions (x and y) of 
//  length "2*rad" and two walls (in z) that are "2*rad" apart
//Input "3" for a cube of length "2*R" that has a spherical 
//  cavity of radius cavRad in the center of it.

//int walls = 1;
//double cavRad;
/*bool periodic_x = false; 
bool periodic_y = false;
bool periodic_z = false;
bool has_x_wall = false;
bool has_y_wall = false;
bool has_z_wall = false;
bool spherical_cavity = false;  //spherical walls on outside of entire volume
bool spherical_solute = true;  //sphere at center of entire volume that is a wall
double lenx = 50; 
double leny = 50; 
double lenz = 50; 
double rad = 50;  //of outer spherical walls
double solRad = 10;  //of inner spherical "solute"*/

//const Vector3d latticex = Vector3d(rad,0,0);

//Vector3d movedj = spheres[j] + Vector3d(rad,0,0);

void run(const double rad, const int N, const long times, const char *filename){
  const double R = 1; 
  Vector3d *spheres = new Vector3d[N];

  // We start with randomly-placed spheres, and then gradually wiggle
  // them around until they are all within bounds and are not
  // overlapping.  We do this by creating an "overlap" value which we
  // constrain to never increase.  Note that this may not work at all
  // for high filling fractions, since we could get stuck in a local
  // minimum.
  
  for(int i=0; i<N; i++){
    spheres[i]=rad*ran3();
  }
  for (int i=0; i<3; i++){
    printf("%f",spheres[1][i]);}
  
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

  FILE *o = fopen(filename, "w");
  fprintf(o,"Radius=%g\n",rad);
  fprintf(o,"Number of Spheres=%d\n", N);
  i = 0;
  long count = 0, workingmoves = 0;

  start = clock();
  num_timed = 0;
  for (long j=0; j<times; j++){
    if (num_timed++ > num_to_time) {
      clock_t now = clock();
      printf("took %g seconds per iteration\n",
             (now - double(start))/CLOCKS_PER_SEC/num_to_time);
      num_timed = 0;
      start = now;
      // after the first timing, just time things once per percent (as
      // often as we print the % complete messages)
      if (times/100 > num_to_time) num_to_time = times/100;
    }
    count++;
    i++;
    // only write out the sphere positions after they've all had a
    // chance to move, to save on file size.
    if (i%N == 0) writeSpheres(spheres, N, o);
    if(j % (times/100)==0){
      printf("%g%% complete...\n",j/(times*1.0)*100);
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
  printf("Total number of attempted moves = %ld\n",count);
  printf("Total number of successful moves = %ld\n",workingmoves);
  printf("Acceptance rate = %g\n", workingmoves/double(count));
  fclose(o);
  delete[] spheres;
}


//const Vector3d latticex = Vector3d(rad,0,0);

//Vector3d movedj = spheres[j] + Vector3d(rad,0,0);

bool periodic_x = false; // will go from -lenx/2 to +lenx/2 
bool periodic_y = false;
bool periodic_z = false;
bool has_x_wall = false;
bool has_y_wall = false;
bool has_z_wall = false;
bool spherical_cavity = false;  //spherical walls on outside of entire volume
bool spherical_solute = true;  //sphere at center of entire volume that is a wall
double lenx = 50; 
double leny = 50; 
double lenz = 50; 
double rad = 50;  //of outer spherical walls
double solRad = 10;  //of inner spherical "solute"
const Vector3d latx = Vector3d(lenx,0,0);
const Vector3d laty = Vector3d(0,leny,0);
const Vector3d latz = Vector3d(0,0,lenz);

double countOverLaps(Vector3d *spheres, int n, double R, double rad){
  double num = 0;
  for(int j = 0; j<n; j++){
    for(int i = j+1; i < n; i++){
        if(distance(spheres[i],spheres[j])<2*R){
          num+=2*R-distance(spheres[i],spheres[j]);
	}
    }
    if (spherical_cavity){
      if(distance(spheres[j],Vector3d(0,0,0))>rad){
        num += distance(spheres[j],Vector3d(0,0,0))-rad;
      }
    }
    if (spherical_solute){
      if(distance(spheres[j],Vector3d(0,0,0))<solRad){
        num -= distance(spheres[j],Vector3d(0,0,0))-solRad;
      }
    }
     //do I need an num += for going outside square walls if move won't allow them to do that?
    if (has_x_wall){
      if (spheres[j][0] > lenx/2 ){
	num += spheres[j][0]-(lenx/2);
      } else if (spheres[j][0] < -lenx/2){
	num -= spheres[j][0] + (lenx/2);
      }
    }
    if (has_y_wall){
      if (spheres[j][1] > leny/2 ){
	num += spheres[j][1]-(leny/2);
      } else if (spheres[j][1] < -leny/2){
	num -= spheres[j][1] + (leny/2);
      }
    }
    if (has_z_wall){
      if (spheres[j][2] > lenz/2 ){
	num += spheres[j][2]-(lenz/2);
      } else if (spheres[j][2] < -lenz/2){
	num -= spheres[j][2] + (lenz/2);
      }
    }
    /*for (int i = j+1; i<n; i++){
      //want to test to make sure that the distance function returns the absolute value (whatever "norm" means)
      if (periodic_x){
        if (distance(spheres[j],(spheres[i]+latx)) < 2*R){
          num += 2*R - distance(spheres[j],(spheres[i]+latx));
	} else if (distance(spheres[j],(spheres[i]-latx)) < 2*R) {
          num += 2*R - distance(spheres[j],(spheres[i]-latx));
	}
      }
      if (periodic_y){
        if (distance(spheres[j],(spheres[i]+laty)) < 2*R){
          num += 2*R - distance(spheres[j],(spheres[i]+laty));
	} else if(distance(spheres[j],(spheres[i]-laty)) < 2*R) {
	  num += 2*R - distance(spheres[j],(spheres[i]-laty));
	}
      }
      if (periodic_z){
        if (distance(spheres[j],(spheres[i]+latz)) < 2*R){
          num += 2*R - distance(spheres[j],(spheres[i]+latz));
	} else if(distance(spheres[j],(spheres[i]-latz)) < 2*R) {
	  num += 2*R - distance(spheres[j],(spheres[i]-latz));
	}
	}*/
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

/*bool periodic_x = false; // will go from -lenx/2 to +lenx/2 
bool periodic_y = false;
bool periodic_z = false;
bool has_x_wall = false;
bool has_y_wall = false;
bool has_z_wall = false;
bool spherical_cavity = false;  //spherical walls on outside of entire volume
bool spherical_solute = true;  //sphere at center of entire volume that is a wall
double lenx = 50; 
double leny = 50; 
double lenz = 50; 
double rad = 50;  //of outer spherical walls
double solRad = 10;  //of inner spherical "solute"
const Vector3d latx = Vector3d(lenx,0,0);
const Vector3d laty = Vector3d(0,leny,0);
const Vector3d latz = Vector3d(0,0,lenz);*/ 

bool overlap(Vector3d *spheres, Vector3d v, int n, double R, double rad, int s){
  for(int i = 0; i < n; i++){
      if(i==s){
        continue;
      }
      if(distance(spheres[i],v)<2*R){
        return true;
      }
      if (periodic_x){
        if (distance(v,spheres[i]+latx) < 2*R){
          return true;
	} else if (distance(v,spheres[i]-latx) < 2*R) {
          return true;
	}
      }
      if (periodic_y){
        if (distance(v,spheres[i]+laty) < 2*R){
          return true;
	} else if (distance(v,spheres[i]-laty) < 2*R) {
          return true;
	}
      }
      if (periodic_z){
        if (distance(v,spheres[i]+latz) < 2*R){
          return true;
	} else if (distance(v,spheres[i]-latz) < 2*R) {
          return true;
	}
      }
  }
    if (spherical_cavity) {
      if (distance(v,Vector3d(0,0,0)) > rad){
        return true;
      }
    }
  if (spherical_solute) {
    if (distance(v,Vector3d(0,0,0)) < solRad){
	return true;
      }
  }
  if (has_x_wall){
    if (v[0] > lenx/2 || v[0] < -lenx/2) {
      return true;
    }
  }
  if (has_y_wall){
    if (v[1] > leny/2 || v[1] < -leny/2) {
      return true;
    }
  }
  if (has_z_wall){
    if (v[2] > lenz/2 || v[2] < -lenz/2) {
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

/*bool periodic_x = false; // will go from -lenx/2 to +lenx/2 
bool periodic_y = false;
bool periodic_z = false;
bool has_x_wall = false;
bool has_y_wall = false;
bool has_z_wall = false;
bool spherical_cavity = false;  //spherical walls on outside of entire volume
bool spherical_solute = true;  //sphere at center of entire volume that is a wall
double lenx = 50; 
double leny = 50; 
double lenz = 50; 
double rad = 50;  //of outer spherical walls
double solRad = 10;  //of inner spherical "solute"
const Vector3d latx = Vector3d(lenx,0,0);
const Vector3d laty = Vector3d(0,leny,0);
const Vector3d latz = Vector3d(0,0,lenz);*/ 


Vector3d move(Vector3d v){
  double scale = .5;
  Vector3d newv = v+scale*ran3();  
  if (periodic_x){
    if (newv[0] > lenx/2){
      newv[0] -= lenx;
    } else if (newv[0] < -lenx/2){
      newv[0] += lenx;
    }
  }
  if (periodic_y){
    if (newv[1] > leny/2){
      newv[1] -= leny;
    } else if (newv[1] < -leny/2){
      newv[1] += leny;
    }
  }
  if (periodic_z){
    if (newv[2] > lenz/2){
      newv[2] -= lenz;
    } else if (newv[2] < -lenz/2){
      newv[2] += lenz;
    }
  }
  return newv;
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
