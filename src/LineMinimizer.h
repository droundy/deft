// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"

typedef Minimizer *(*LineMinimizer)(Functional f, VectorXd *data, const VectorXd &direction,
                                    double gradDotDirection, double *step);

Minimizer *QuadraticLineMinimizer(Functional f, VectorXd *data, const VectorXd &direction,
                                  double gradDotDirection, double *step);
