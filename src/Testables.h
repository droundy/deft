// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"

Functional eta_effective(Functional eta, double lambdainput);
Functional gSW(double temp, double R, double epsdis0, double lambda);
Functional da1_dlam(double radius, double epsdis, double lambdainput);
