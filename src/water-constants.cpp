// Deft is a density functional package developed by the research
// group of Professor David Roundy
//
// Copyright 2010 The Deft Authors
//
// Deft is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with deft; if not, write to the Free Software Foundation,
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// Please see the file AUTHORS for a list of authors.

#include "equation-of-state.h"

static const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin
static const double angstrom = 1.8897261; // An angstrom in atomic units

LiquidProperties hughes_water_prop = {
  //2.93, // hard sphere radius of water in bohr
  3.03420*angstrom/2, // hard sphere radius of water in bohr from Clark et al, table 6
  4.9388942e-3, // density of liquid water
  0.0017, // approximate critical density of water
  1.14e-7, // vapor density of water
  298.15*kB, // room temperature in Hartree
  250.0*kB, // epsilon_dispersion from Clark et al, table 6
  1.78890, // lambda_dispersion from Clark et al, table 6
  0.353, // length-scaling parameter that DJR added to fix surface tension
  1400.00*kB, // epsilonAB from Clark et al, table 6

  // The following two are almost the same, but differ ever so
  // slightly.  Using the decimal rounded value from the paper gives
  // ever so small differences from the value actually used (which is
  // a round decimal relative to \sigma^3.
  0.0381876*(3.03420*angstrom*3.03420*angstrom*3.03420*angstrom), // kappaAB from Clark et al code
  //1.06673*angstrom*angstrom*angstrom, // kappaAB from Clark et al, table 6
};

LiquidProperties new_water_prop = {
  //2.93, // hard sphere radius of water in bohr
  3.03420*angstrom/2, // hard sphere radius of water in bohr from Clark et al, table 6
  4.9388942e-3, // density of liquid water
  0.0017, // approximate critical density of water
  1.14e-7, // vapor density of water
  298.15*kB, // room temperature in Hartree
  250.0*kB, // epsilon_dispersion from Clark et al, table 6
  1.78890, // lambda_dispersion from Clark et al, table 6
  0.454, // length-scaling parameter that DJR added to fix surface tension
  1400.00*kB, // epsilonAB from Clark et al, table 6

  // The following two are almost the same, but differ ever so
  // slightly.  Using the decimal rounded value from the paper gives
  // ever so small differences from the value actually used (which is
  // a round decimal relative to \sigma^3.
  0.0381876*(3.03420*angstrom*3.03420*angstrom*3.03420*angstrom), // kappaAB from Clark et al code
  //1.06673*angstrom*angstrom*angstrom, // kappaAB from Clark et al, table 6
};
