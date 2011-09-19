#pragma once

#include <Eigen/Core>
#include "MersenneTwister.h"
USING_PART_OF_NAMESPACE_EIGEN

double ran();
Vector3d ran3();
double distance(Vector3d v1, Vector3d v2);
void writeSpheres(Vector3d *spheres, int n, FILE * o);
Vector3d move(Vector3d v, double x, double y, double z);
void run(const double rad, const int N, const int times);
bool overlap(Vector3d *spheres, Vector3d v, int n, double R, int s);
bool touch(Vector3d *spheres, Vector3d v, int n, double R, double delta, int s);
Vector3d move(Vector3d v, double R);
Vector3d move(Vector3d v);
bool overlap(Vector3d *spheres, Vector3d v, int n, double R, double rad, int s);
bool overlap(Vector3d *spheres, Vector3d v, int n, double R, int s, double x, double y, double z);
int shell(Vector3d v, int div, double R);
double countOverLaps(Vector3d *spheres, int n, double R, double rad);
