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
#include "equation-of-state.h"
#include "Minimizer.h"
#include "Functionals.h"
#include <stdio.h>

long Eigen::djr_memused = 0;
long Eigen::djr_mempeak = 0;

long peak_memory() {
  return Eigen::djr_mempeak;
}

long current_memory() {
  return Eigen::djr_memused;
}

void reset_peak_memory() {
  Eigen::djr_mempeak = Eigen::djr_memused;
}
