#pragma once

#include <Eigen/Core>
#include "MersenneTwister.h"
USING_PART_OF_NAMESPACE_EIGEN

double ran();
Vector3d ran3();
void writeSpheres(Vector3d *spheres, int n, FILE * o);
Vector3d move(Vector3d v, double x, double y, double z);
void run(const double rad, const int N, const long times, const char *filename);
bool overlap(Vector3d *spheres, Vector3d v, int n, double R, int s);
bool touch(Vector3d w, Vector3d v, double oShell);
Vector3d move(Vector3d v, double R);
Vector3d move(Vector3d v);
int shell(Vector3d v, int div, double *radius, double *sections);
double countOverLaps(Vector3d *spheres, int n, double R);

inline double distance(Vector3d v1, Vector3d v2){
  return (v1 - v2).norm();
}
