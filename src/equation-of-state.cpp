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
#include "Functionals.h"
#include <stdio.h>

double find_minimum(Functional f, double nmin, double nmax) {
  double nbest = nmin;
  double ebest = f(nmin);
  //printf("Limits are %g and %g\n", nmin, nmax);
  const double dn = (nmax - nmin)*1e-2;
  for (double n = nmin; n<=nmax; n += dn) {
    double en = f(n);
    // printf("Considering %g with energy %g\n", n, en);
    if (en < ebest) {
      ebest = en;
      nbest = n;
    }
  }
  //printf("best Veff is %g\n", nbest);
  double nlo = nbest - dn, nhi = nbest + dn;
  while ((nhi - nlo)/fabs(nbest) > 1e-15) {
    if (nbest < 0.5*(nhi+nlo)) {
      double ntry = 0.3*nlo + 0.7*nhi;
      double etry = f(ntry);
      if (etry < ebest) {
        nlo = nbest;
        nbest = ntry;
        ebest = etry;
      } else {
        nhi = ntry;
      }
    } else {
      double ntry = 0.7*nlo + 0.3*nhi;
      double etry = f(ntry);
      if (etry < ebest) {
        nhi = nbest;
        nbest = ntry;
        ebest = etry;
      } else {
        nlo = ntry;
      }
    }
    //printf("improved Veff is %g\n", nbest);
  }
  return nbest;
}

double pressure(Functional f, double kT, double density) {
  double V = -kT*log(density);
  //printf("density is %g\n", density);
  //printf("f(V) is %g\n", f(V));
  //printf("-f.derive(V)*kT is %g\n", -f.derive(V)*kT);
  return -f.derive(V)*kT - f(V);
}

double pressure_to_density(Functional f, double kT, double p, double nmin, double nmax) {
  while (nmax/nmin > 1 + 1e-14) {
    double ntry = sqrt(nmax*nmin);
    double ptry = pressure(f, kT, ntry);
    if (ptry < p) nmin = ntry;
    else nmax = ntry;
  }
  return sqrt(nmax*nmin);
}

double find_density(Functional f, double kT, double nmin, double nmax) {
  double Vmax = -kT*log(nmin);
  double Vmin = -kT*log(nmax);
  double V = find_minimum(f, Vmin, Vmax);
  return EffectivePotentialToDensity(kT)(V);
}

double find_chemical_potential(Functional f, double kT, double n) {
  const double V = -kT*log(n);
  return -f.derive(V)*kT/n;
}

void other_equation_of_state(FILE *o, Functional f, double kT, double nmin, double nmax) {
  Functional n = EffectivePotentialToDensity(kT);
  const double mumin = find_chemical_potential(f, kT, nmin);
  const double mumax = find_chemical_potential(f, kT, nmax);
  const double dmu = (mumax - mumin)/2000;
  printf("mumin is %g and mumax is %g\n", mumin, mumax);
  for (double mu=mumin; mu<=mumax; mu += dmu) {
    printf("Working on mu = %g\n", mu);
    Functional fmu = f + ChemicalPotential(mu)(n);
    double n = find_density(fmu, kT, nmin, nmax);
    printf("Got density = %g\n", n);
    double V = -kT*log(n);
    double p = pressure(fmu, kT, n);
    //printf("Got ngoal of %g but mu is %g\n", ngoal, mu);
    double e = f(V);
    fprintf(o, "%g\t%g\t%g\n", n, p, e);
  }
}

void equation_of_state(FILE *o, Functional f, double kT, double nmin, double nmax) {
  const double factor = 1.04;
  for (double ngoal=nmin; ngoal<nmax; ngoal *= factor) {
    double V = -kT*log(ngoal);
    double p = pressure(f, kT, ngoal);
    double der = -f.derive(V)*kT/ngoal;
    //printf("Got ngoal of %g but mu is %g\n", ngoal, mu);
    double e = f(V);
    fprintf(o, "%g\t%g\t%g\t%g\n", ngoal, p, e, der);
  }
}
