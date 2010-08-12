// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"

typedef Minimizer *(*LineMinimizer)(Functional f, const GridDescription &gd, VectorXd *data,
                                    const VectorXd &direction, double gradDotDirection, double *step);

Minimizer *QuadraticLineMinimizer(Functional f, const GridDescription &gd, VectorXd *data,
                                  const VectorXd &direction, double gradDotDirection, double *step);
