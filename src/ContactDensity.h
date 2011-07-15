// -*- mode: C++; -*-

#pragma once

#include "Functionals.h"

Functional dWBNT_dn0(double radius);
Functional dWBNT_dn1(double radius);
Functional dWBNT_dn2(double radius);
Functional dWBNT_dn3(double radius);

Functional dWBNT_dnV1(double radius);
Functional dWBNT_dnV2(double radius);

Functional ContactDensitySimplest(double radius);
Functional ContactDensitySphere(double radius);

Functional FuWuContactDensity(double radius);
