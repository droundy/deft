#pragma once

#include <Eigen/Core>

USING_PART_OF_NAMESPACE_EIGEN

double ran();
Vector3d ran3();
double distance(Vector3d v);
void writeSpheres(Vector3d *spheres, int n, FILE * o);
