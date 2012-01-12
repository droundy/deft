#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"

//Input "1" for large spherical cavity of radius "rad"
//Input "2" for two dimensions boundary conditions (x and y) of 
//  length "2*rad" and two walls (in z) that are "2*rad" apart
//Input "3" for a cube of length "2*R" that has a spherical 
//  cavity of radius cavRad in the center of it.

int walls = 1;
double cavRad;


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
    spheres[i%N]=move(spheres[i%N],rad);
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
    Vector3d temp = move(spheres[i%N], rad);
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


double countOverLaps(Vector3d *spheres, int n, double R, double rad){
  double num = 0;
  if (walls == 1){
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
  } else if (walls ==2){
    //the following ignores two spheres in opposite corners of large cube
    for(int j = 0; j<n; j++){
      if (spheres[j][2] > rad){
        num += (spheres[j][2] - rad);
      } else if (spheres[j][2] < (-1)*rad){
        num -= (spheres[j][2] + rad);
      }
      for (int i = j+1; i<n; i++){
        if(distance(spheres[j],Vector3d(0,0,0))>rad){
          num+=distance(spheres[j],Vector3d(0,0,0))-rad;
        }
        if (rad-spheres[j][0] < R || spheres[j][0]-rad < R ||rad-spheres[j][1] < R || spheres[j][1]-rad < R){  
          if ( (spheres[j][0]-spheres[i][0]+2*rad)*(spheres[j][0]-spheres[i][0]+2*rad) 
	       + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
	    num -= (spheres[j][0]-spheres[i][0]+2*rad)*(spheres[j][0]-spheres[i][0]+2*rad) 
		     + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
		     + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) - 2*rad*2*rad;
	  } 
          if ( (spheres[j][0]-spheres[i][0]-2*rad)*(spheres[j][0]-spheres[i][0]-2*rad) 
	       + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1])
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
            num -= (spheres[j][0]-spheres[i][0]+2*rad)*(spheres[j][0]-spheres[i][0]+2*rad) 
		     + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
		     + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) - 2*rad*2*rad;
          }
          if ( (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
	       + (spheres[j][1]-spheres[i][1]+2*rad)*(spheres[j][1]-spheres[i][1]+2*rad) 
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
            num -= (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
		     + (spheres[j][1]-spheres[i][1]+2*rad)*(spheres[j][1]-spheres[i][1]+2*rad) 
		     + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) - 2*rad*2*rad;
	  }
          if ( (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
	       + (spheres[j][1]-spheres[i][1]-2*rad)*(spheres[j][1]-spheres[i][1]-2*rad) 
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
            num -= (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
		     + (spheres[j][1]-spheres[i][1]-2*rad)*(spheres[j][1]-spheres[i][1]-2*rad) 
		     + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) - 2*rad*2*rad;
	  }
	}
      }
    }
  } else if (walls == 3){
    for(int j = 0; j<n; j++){
      if(distance(spheres[j],Vector3d(0,0,0)) < cavRad){
        num -= distance(spheres[j],Vector3d(0,0,0))-cavRad;
      }
      for(int i = j+1; i < n; i++){
        if(distance(spheres[i],spheres[j])<2*R){
          num+=2*R-distance(spheres[i],spheres[j]);
	}
        if (rad-spheres[j][0] < R || spheres[j][0]-rad < R || rad-spheres[j][1] < R || 
            spheres[j][1]-rad < R ||rad-spheres[j][2] < R || spheres[j][2]-rad < R ){  
          if ( (spheres[j][0]-spheres[i][0]+2*rad)*(spheres[j][0]-spheres[i][0]+2*rad) 
	       + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
	    num -= (spheres[j][0]-spheres[i][0]+2*rad)*(spheres[j][0]-spheres[i][0]+2*rad) 
		     + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
		     + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) - 2*rad*2*rad;
	  } 
          if ( (spheres[j][0]-spheres[i][0]-2*rad)*(spheres[j][0]-spheres[i][0]-2*rad) 
	       + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1])
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
            num -= (spheres[j][0]-spheres[i][0]+2*rad)*(spheres[j][0]-spheres[i][0]+2*rad) 
		     + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
		     + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) - 2*rad*2*rad;
          }
          if ( (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
	       + (spheres[j][1]-spheres[i][1]+2*rad)*(spheres[j][1]-spheres[i][1]+2*rad) 
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
            num -= (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
		     + (spheres[j][1]-spheres[i][1]+2*rad)*(spheres[j][1]-spheres[i][1]+2*rad) 
		     + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) - 2*rad*2*rad;
	  }
          if ( (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
	       + (spheres[j][1]-spheres[i][1]-2*rad)*(spheres[j][1]-spheres[i][1]-2*rad) 
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
            num -= (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
		     + (spheres[j][1]-spheres[i][1]-2*rad)*(spheres[j][1]-spheres[i][1]-2*rad) 
		     + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) - 2*rad*2*rad;
	  }
          if ( (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
	       + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
	       + (spheres[j][2]-spheres[i][2]+2*rad)*(spheres[j][2]-spheres[i][2]+2*rad) < 2*rad*2*rad){
            num -= (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
		     + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
		     + (spheres[j][2]-spheres[i][2]+2*rad)*(spheres[j][2]-spheres[i][2]+2*rad) - 2*rad*2*rad;
	  }
          if ( (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
	       + (spheres[j][1]-spheres[i][1]-2*rad)*(spheres[j][1]-spheres[i][1]-2*rad) 
	       + (spheres[j][2]-spheres[i][2])*(spheres[j][2]-spheres[i][2]) < 2*rad*2*rad){
            num -= (spheres[j][0]-spheres[i][0])*(spheres[j][0]-spheres[i][0]) 
		     + (spheres[j][1]-spheres[i][1])*(spheres[j][1]-spheres[i][1]) 
		     + (spheres[j][2]-spheres[i][2]-2*rad)*(spheres[j][2]-spheres[i][2]-2*rad) - 2*rad*2*rad;
	  }
	}
      }
    }
  }
 return num;
}

bool overlap(Vector3d *spheres, Vector3d v, int n, double R, double rad, int s){
  for(int i = 0; i < n; i++){
      if(i==s){
        continue;
      }
      if(distance(spheres[i],v)<2*R){
        return true;
      }
  }
  for(int i = 0; i < n; i++){
    if (walls == 1){
      if(distance(v,Vector3d(0,0,0))>rad){
          return true;
      }
    } else if (walls == 2){
    if (v[2] > rad || v[2] < (-1)*rad){
      return true;
    }
    if ( (v[0]-spheres[i][0]+2*rad)*(v[0]-spheres[i][0]+2*rad) 
         + (v[1]-spheres[i][1])*(v[1]-spheres[i][1]) 
         + (v[2]-spheres[i][2])*(v[2]-spheres[i][2]) < 2*rad*2*rad ||
         (v[0]-spheres[i][0]-2*rad)*(v[0]-spheres[i][0]-2*rad) 
         + (v[1]-spheres[i][1])*(v[1]-spheres[i][1])
         + (v[2]-spheres[i][2])*(v[2]-spheres[i][2]) < 2*rad*2*rad ||
         (v[0]-spheres[i][0])*(v[0]-spheres[i][0]) 
         + (v[1]-spheres[i][1]+2*rad)*(v[1]-spheres[i][1]+2*rad) 
         + (v[2]-spheres[i][2])*(v[2]-spheres[i][2]) < 2*rad*2*rad ||
         (v[0]-spheres[i][0])*(v[0]-spheres[i][0]) 
         + (v[1]-spheres[i][1]-2*rad)*(v[1]-spheres[i][1]-2*rad) 
         + (v[2]-spheres[i][2])*(v[2]-spheres[i][2]) < 2*rad*2*rad ){
      return true;}
    } else if (walls == 3){
      if (distance(v,Vector3d(0,0,0)) < cavRad){
      return true;
    }
      if ( (v[0]-spheres[i][0]+2*rad)*(v[0]-spheres[i][0]+2*rad) 
           + (v[1]-spheres[i][1])*(v[1]-spheres[i][1]) 
           + (v[2]-spheres[i][2])*(v[2]-spheres[i][2]) < 2*rad*2*rad ||
         (v[0]-spheres[i][0]-2*rad)*(v[0]-spheres[i][0]-2*rad) 
           + (v[1]-spheres[i][1])*(v[1]-spheres[i][1])
           + (v[2]-spheres[i][2])*(v[2]-spheres[i][2]) < 2*rad*2*rad ||
         (v[0]-spheres[i][0])*(v[0]-spheres[i][0]) 
           + (v[1]-spheres[i][1]+2*rad)*(v[1]-spheres[i][1]+2*rad) 
           + (v[2]-spheres[i][2])*(v[2]-spheres[i][2]) < 2*rad*2*rad ||
         (v[0]-spheres[i][0])*(v[0]-spheres[i][0]) 
           + (v[1]-spheres[i][1]-2*rad)*(v[1]-spheres[i][1]-2*rad) 
	   + (v[2]-spheres[i][2])*(v[2]-spheres[i][2]) < 2*rad*2*rad ||
         (v[0]-spheres[i][0])*(v[0]-spheres[i][0]) 
	   + (v[1]-spheres[i][1])*(v[1]-spheres[i][1]) 
           + (v[2]-spheres[i][2]+2*rad)*(v[2]-spheres[i][2]+2*rad) < 2*rad*2*rad ||
         (v[0]-spheres[i][0])*(v[0]-spheres[i][0]) 
	   + (v[1]-spheres[i][1])*(v[1]-spheres[i][1]) 
	   + (v[2]-spheres[i][2]-2*rad)*(v[2]-spheres[i][2]-2*rad) < 2*rad*2*rad){
      return true;}
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


Vector3d move(Vector3d v, double rad){
  double scale = .5;
  if (walls == 1){
    return v+scale*ran3();
  } else if (walls == 2){
    Vector3d newv = v +scale*ran3();
    if (newv[0] > rad){
      newv[0] -= 2*rad;
    }  else if (newv[0] < (-1)*rad){
      newv[0] += 2*rad;}
    if (newv[1] > rad){
      newv[1] -= 2*rad;
    }  else if (newv[1] < (-1)*rad){
      newv[1] += 2*rad;}
    return newv;
  } else if (walls == 3){
    Vector3d newv = v +scale*ran3();
    if (newv[0] > rad){
      newv[0] -= 2*rad;
    }  else if (newv[0] < (-1)*rad){
      newv[0] += 2*rad;}
    if (newv[1] > rad){
      newv[1] -= 2*rad;
    }  else if (newv[1] < (-1)*rad){
      newv[1] += 2*rad;}
    if (newv[2] > rad){
      newv[2] -= 2*rad;
    }  else if (newv[2] < (-1)*rad){
      newv[2] += 2*rad;}
    return newv;
  }
  printf("You need to pick a Boundary type!");
  fflush(stdout);
  exit (1);
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
