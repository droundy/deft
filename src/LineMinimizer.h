// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"

typedef Minimizer (*LineMinimizer)(FieldFunctional f, const GridDescription &gd, VectorXd *data,
                                    const VectorXd &direction, double gradDotDirection, double *step);

Minimizer QuadraticLineMinimizer(FieldFunctional f, const GridDescription &gd, VectorXd *data,
                                 const VectorXd &direction, double gradDotDirection, double *step);

Minimizer SteepestDescent(FieldFunctional f, const GridDescription &gdin, VectorXd *data,
                          LineMinimizer lm, double stepsize = 0.1);
Minimizer PreconditionedSteepestDescent(FieldFunctional f, const GridDescription &gdin, VectorXd *data,
                                        LineMinimizer lm, double stepsize = 0.1);

Minimizer ConjugateGradient(FieldFunctional f, const GridDescription &gdin, VectorXd *data,
                            LineMinimizer lm, double stepsize = 0.1);
Minimizer PreconditionedConjugateGradient(FieldFunctional f, const GridDescription &gdin, VectorXd *data,
                                          LineMinimizer lm, double stepsize = 0.1);
