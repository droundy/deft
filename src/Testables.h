// -*- mode: C++; -*-

#pragma once

#include "Grid.h"
#include "Functional.h"

Functional eta_effective(Functional eta, double lambdainput);
Functional gSW(double R, double epsdis0, double lambda, double lscale);
Functional da1_dlam(double radius, double epsdis, double lambdainput, double lscale);
Functional da1_deta(double radius, double epsdis, double lambdainput, double lscale);
