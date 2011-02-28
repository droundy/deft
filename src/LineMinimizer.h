// -*- mode: C++; -*-

#pragma once

#include "Minimizer.h"

typedef Minimizer (*LineMinimizer)(Functional f, const GridDescription &gd, double kT, VectorXd *data,
                                    const VectorXd &direction, double gradDotDirection, double *step);

Minimizer QuadraticLineMinimizer(Functional f, const GridDescription &gd, double kT, VectorXd *data,
                                 const VectorXd &direction, double gradDotDirection, double *step);

Minimizer SteepestDescent(Functional f, const GridDescription &gdin, double kT, VectorXd *data,
                          LineMinimizer lm, double stepsize = 0.1);
Minimizer PreconditionedSteepestDescent(Functional f, const GridDescription &gdin, double kT, VectorXd *data,
                                        LineMinimizer lm, double stepsize = 0.1);

Minimizer ConjugateGradient(Functional f, const GridDescription &gdin, double kT, VectorXd *data,
                            LineMinimizer lm, double stepsize = 0.1);
Minimizer PreconditionedConjugateGradient(Functional f, const GridDescription &gdin, double kT, VectorXd *data,
                                          LineMinimizer lm, double stepsize = 0.1);
