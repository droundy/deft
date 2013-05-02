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

Functional gSigmaA(double radius);
Functional gSigmaS(double radius);
Functional gSigmaS2(double radius);
Functional gSigmaA2(double radius);
Functional gSigmaS_m2(double radius); // from Haskell
Functional gSigmaA_m2(double radius); // from Haskell
Functional gSigmaAm2(double radius);
Functional gSigmaSm2(double radius);
Functional CorrelationGrossCorrect(double radius);

Functional GrossCorrelation(double radius);
Functional YuWuCorrelation(double radius);
Functional YuWuCorrelation_S(double radius);
