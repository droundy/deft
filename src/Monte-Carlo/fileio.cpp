#include <stdio.h>
#include "monte-carlo.h"

void writeSpheres(Vector3d *spheres, int n, FILE * o){
  for(int i = 0; i < n; i++){
    fprintf(o,"%g\t%g\t%g\t", spheres[i][0], spheres[i][1], spheres[i][2]);
  }
  fprintf(o,"\n");
  fflush(o);
}
