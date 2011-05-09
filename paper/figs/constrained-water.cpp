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

#include <stdio.h>
#include <time.h>
#include "Functionals.h"
#include "equation-of-state.h"
#include "LineMinimizer.h"

#include "utilities.h"

// Here we set up the lattice.
const double zmax = 80;
const double width = 0.0001;
const double cavitysize = 40;

double notinwall(Cartesian r) {
  const double z = r.dot(Cartesian(0,0,1));
  if (fabs(z) > cavitysize/2) {
    return 0;
  } else {
    return 1;
  }
}

int main(int, char **) {
  FILE *o = fopen("paper/figs/constrained-water.dat", "w");

  Functional f = OfEffectivePotential(SaftFluidSlow(water_prop.lengthscale,
                                                    water_prop.epsilonAB, water_prop.kappaAB,
                                                    water_prop.epsilon_dispersion,
                                                    water_prop.lambda_dispersion, water_prop.length_scaling, 0));
  double mu_satp = find_chemical_potential(f, water_prop.kT,
                                           water_prop.liquid_density);
  f = OfEffectivePotential(SaftFluidSlow(water_prop.lengthscale,
                                         water_prop.epsilonAB, water_prop.kappaAB,
                                         water_prop.epsilon_dispersion,
                                         water_prop.lambda_dispersion, water_prop.length_scaling, mu_satp));

  Lattice lat(Cartesian(width,0,0), Cartesian(0,width,0), Cartesian(0,0,zmax));
  GridDescription gd(lat, 0.1);

  Grid potential(gd);
  Grid constraint(gd);
  constraint.Set(notinwall);
  //f = constrain(constraint, f);

  Minimizer min = Precision(0, ConjugateGradient(f, gd, water_prop.kT, &potential,
						 QuadraticLineMinimizer));

  potential = water_prop.liquid_density*constraint
    + water_prop.vapor_density*VectorXd::Ones(gd.NxNyNz);
  //potential = water_prop.liquid_density*VectorXd::Ones(gd.NxNyNz);
  potential = -water_prop.kT*potential.cwise().log();

  const int numiters = 50;
  for (int i=0;i<numiters && min.improve_energy(true);i++) {
    fflush(stdout);
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    density.epsNative1d("paper/figs/1d-constrained-plot.eps",
			Cartesian(0,0,0), Cartesian(0,0,zmax),
			water_prop.liquid_density, 1, "Y axis: , x axis: ");
  }
  min.print_info();

  double energy = min.energy()/width/width;
  printf("Energy is %.15g\n", energy);

  double N = 0;
  {
    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    for (int i=0;i<gd.NxNyNz;i++) N += density[i]*gd.dvolume;
  }
  N = N/width/width;
  printf("N is %.15g\n", N);

    Grid density(gd, EffectivePotentialToDensity()(water_prop.kT, gd, potential));
    density.epsNative1d("paper/figs/1d-constrained-plot.eps", Cartesian(0,0,0), Cartesian(0,0,zmax), water_prop.liquid_density, 1, "Y axis: , x axis: ");
    //potential.epsNative1d("hard-wall-potential.eps", Cartesian(0,0,0), Cartesian(0,0,zmax), 1, 1);

    fclose(o);
}
