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

#include "utilities.h"

static const double kB = 3.16681539628059e-6; // This is Boltzmann's constant in Hartree/Kelvin

LiquidProperties water_prop = {
  2.93, // hard sphere radius of water in bohr
  4.9388942e-3, // density of liquid water
  1.14e-7, // vapor density of water
  298.15*kB, // room temperature in Hartree
  290.1*kB, // epsilon_dispersion
  2859*kB, // epsilonAB
  0.023, // kappaAB
};
