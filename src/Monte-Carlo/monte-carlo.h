#pragma once

#include <Eigen/Core>

USING_PART_OF_NAMESPACE_EIGEN

double ran();
Vector3d ran3();
double distance(Vector3d v1, Vector3d v2);
void writeSpheres(Vector3d *spheres, int n, FILE * o);
Vector3d move(Vector3d v, double x, double y, double z);
void run();
bool overlap(Vector3d *spheres, Vector3d v, int n, double R, int s);
bool touch(Vector3d *spheres, Vector3d v, int n, double R, double delta, int s);
Vector3d move(Vector3d v, double R);
