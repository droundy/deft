// -*- mode: C++; -*-

#pragma once

#include "Functionals.h"

Functional dWBNT_dn0(double radius);
Functional dWBNT_dn1(double radius);
Functional dWBNT_dn2(double radius);
Functional dWBNT_dn3(double radius);

Functional dWBNT_dnV1(double radius);
Functional dWBNT_dnV2(double radius);

Functional dWBm2_dn0(double radius);
Functional dWBm2_dn1(double radius);
Functional dWBm2_dn2(double radius);
Functional dWBm2_dn3(double radius);

Functional dWBm2_dnV1(double radius);
Functional dWBm2_dnV2(double radius);

Functional dAdR_simplest(double radius);

Functional dAdR_S_WBNT(double radius);
Functional dAdR_sphere_over_n_WBNT(double radius);

Functional dAdR_S_WBm2(double radius);
Functional dAdR_sphere_over_n_WBm2(double radius);

Functional ContactDensitySimplest(double radius);
Functional ContactDensitySphere(double radius);
Functional ContactDensity_S(double radius);
Functional ContactDensitySphereWBm2(double radius);
Functional ContactDensity_S_WBm2(double radius);

Functional GrossContactDensity(double radius);
Functional FuWuContactDensity(double radius);
Functional FuWuContactDensity_S(double radius);
