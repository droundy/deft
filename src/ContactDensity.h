// -*- mode: C++; -*-

#pragma once

#include "Functionals.h"

Functional n0(double radius);

Functional dWB_dn0(double radius);
Functional dWB_dn1(double radius);
Functional dWB_dn2(double radius);
Functional dWB_dn3(double radius);

Functional dWB_dn1v_over_n2v(double radius);
Functional dWB_dn2v_over_n2v(double radius);

Functional phi2(const Functional &n3);
Functional phi3(const Functional &n3);

Functional dWBm2_dn0(double radius);
Functional dWBm2_dn1(double radius);
Functional dWBm2_dn2(double radius);
Functional dWBm2_dn3(double radius);

Functional dWBm2_dn1v_over_n2v(double radius);
Functional dWBm2_dn2v_over_n2v(double radius);

Functional dAdR_S(double radius);
Functional dAdR_A_over_n(double radius);

Functional dAdR_S_WBm2(double radius);
Functional dAdR_A_over_n_WBm2(double radius);

Functional Correlation_A(double radius);
Functional Correlation_S(double radius);
Functional Correlation_A_WBm2(double radius);
Functional Correlation_S_WBm2(double radius);

Functional GrossCorrelation(double radius);
Functional YuWuCorrelation(double radius);
Functional YuWuCorrelation_S(double radius);
