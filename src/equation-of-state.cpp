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

static inline double sqr(double x) {
  return x*x;
}

double find_minimum_slow(Functional f, double kT, double xmin, double xmax) {
  double dmin = f.derive(kT, xmin);
  double dmax = f.derive(kT, xmax);
  if (dmax < 0) {
    printf("Oh no, xmax has negative slope!\n");
    exit(1);
  }
  if (dmin > 0) {
    printf("Oh no, xmin has positive slope!\n");
    exit(1);
  }
  const double fraccuracy = 1e-15;
  do {
    // First let's do a bisection (to guarantee we converge eventually).
    double xtry = 0.5*(xmax + xmin);
    double dtry = f.derive(kT, xtry);
    if (dtry > 0) {
      dmax = dtry;
      xmax = xtry;
    } else if (dtry == 0) {
      return xtry;
    } else {
      xmin = xtry;
      dmin = dtry;
    }
    // Now let's do a secant method (so maybe we'll converge really quickly)...
    xtry = (xmin*dmax - xmax*dmin)/(dmax - dmin);
    dtry = f.derive(kT, xtry);
    if (dtry > 0) {
      dmax = dtry;
      xmax = xtry;
    } else if (dtry == 0) {
      return xtry;
    } else {
      xmin = xtry;
      dmin = dtry;
    }
    // Now let's use a secant approach to "bracket" the minimum hopefully.
    xtry = (xmin*dmax - xmax*dmin)/(dmax - dmin);
    if (xtry - xmin > xmax - xtry) {
      xtry -= xmax - xtry;
    } else {
      xtry += xtry - xmin;
    }
    dtry = f.derive(kT, xtry);
    if (dtry > 0) {
      dmax = dtry;
      xmax = xtry;
    } else if (dtry == 0) {
      return xtry;
    } else {
      xmin = xtry;
      dmin = dtry;
    }
  } while (fabs((xmax-xmin)/xmax) > fraccuracy);
  return 0.5*(xmax+xmin);
}

double find_minimum(Functional f, double kT, double nmin, double nmax) {
  double nbest = nmin;
  double ebest = f(kT, nmin);
  //printf("Limits are %g and %g\n", nmin, nmax);
  const double dn = (nmax - nmin)*1e-2;
  for (double n = nmin; n<=nmax; n += dn) {
    double en = f(kT, n);
    // printf("Considering %g with energy %g\n", n, en);
    if (en < ebest) {
      ebest = en;
      nbest = n;
    }
  }
  //printf("best Veff is %g\n", nbest);
  double nlo = nbest - dn, nhi = nbest + dn;
  double elo = f(kT, nlo);
  double ehi = f(kT, nhi);
  const double fraccuracy = 1e-15;
  while ((nhi - nlo)/fabs(nbest) > fraccuracy) {
    // First we'll do a golden-section search...
    if (nbest < 0.5*(nhi+nlo)) {
      double ntry = 0.38*nlo + 0.62*nhi;
      double etry = f(kT, ntry);
      if (etry < ebest) {
        nlo = nbest;
        elo = ebest;
        nbest = ntry;
        ebest = etry;
      } else {
        nhi = ntry;
        ehi = etry;
      }
    } else {
      double ntry = 0.62*nlo + 0.38*nhi;
      double etry = f(kT, ntry);
      if (etry < ebest) {
        nhi = nbest;
        ehi = ebest;
        nbest = ntry;
        ebest = etry;
      } else {
        nlo = ntry;
        elo = etry;
      }
    }

    // Now let's try something a bit more bold...
    
    double ntry = nbest - 0.5*(sqr(nbest-nlo)*(ebest-elo) - sqr(nbest-nhi)*(ebest-ehi))/
      ((nbest-nlo)*(ebest-elo) - (nbest-nhi)*(ebest-ehi));
    if (ntry > nlo && ntry < nhi && ntry != nbest) {
      double etry = f(kT, ntry);
      if (ntry < nbest) {
        if (etry < ebest) {
          nhi = nbest;
          ehi = ebest;
          nbest = ntry;
          ebest = etry;
        } else {
          nlo = ntry;
          elo = etry;
        }
      } else {
        if (etry < ebest) {
          nlo = nbest;
          elo = ebest;
          nbest = ntry;
          ebest = etry;
        } else {
          nhi = ntry;
          ehi = etry;
        }
      }
    }
    
    //printf("improved Veff is %g +/- %g\n", nbest, (nhi-nlo));
  }
  return nbest;
}

// p = n df/dn - f = -kT df/dVeff - f
double pressure(Functional f, double kT, double density) {
  double V = -kT*log(density);
  //printf("density is %g\n", density);
  //printf("f(V) is %g\n", f(V));
  //printf("-f.derive(V)*kT is %g\n", -f.derive(V)*kT);
  return -f.derive(kT, V)*kT - f(kT, V);
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
  double V = find_minimum(f, kT, Vmin, Vmax);
  return EffectivePotentialToDensity()(kT, V);
}

double find_chemical_potential(Functional f, double kT, double n) {
  const double V = -kT*log(n);
  return f.derive(kT, V)*kT/n;
}

double chemical_potential_to_density(Functional f, double kT, double mu,
                                     double nmin, double nmax) {
  Functional n = EffectivePotentialToDensity();
  Functional fmu = f + ChemicalPotential(mu)(n);
  return find_density(fmu, kT, nmin, nmax);
}

double coexisting_vapor_density(Functional f, double kT, double liquid_density,
                                double vapor_liquid_ratio) {
  const double mu = find_chemical_potential(f, kT, liquid_density);
  double nmax = liquid_density;
  // Here I assume that the vapor density is at least a hundredth of a
  // percent smaller than the liquid density.
  Functional n = EffectivePotentialToDensity();
  Functional fmu = f + ChemicalPotential(mu)(n);

  // The following is a very rough guess as to the density of the
  // vapor based on the idea that it is an ideal gas with pressure
  // equal to the pressure of the liquid.  We look at twice this
  // density and see if this is higher than the liquid density.
  double ng_simple = pressure(f, kT, liquid_density)/kT;
  if (ng_simple <= 0) ng_simple = liquid_density*1e-6;
  if (-find_chemical_potential(fmu, kT, ng_simple*2) > 0 &&
      ng_simple < 0.5*liquid_density && false) {
    nmax = ng_simple*2;
  } else {
    nmax = liquid_density;
    const double root_vapor_liquid_ratio = sqrt(vapor_liquid_ratio);
    double deriv;
    do {
      // This is a cautious approach which should work up to the point
      // where the liquid density and vapor density are within something
      // like vapor_liquid_ratio of each other.  Thus as you approach
      // the critical point you need to adjust vapor_liquid_ratio.  We
      // use root_vapor_liquid_ratio which is closer to 1 so that we've
      // still got some wiggle room.
      nmax *= root_vapor_liquid_ratio;
      deriv = -find_chemical_potential(fmu, kT, nmax);
    } while (deriv < 0 && nmax > 0.0001);
  }

  //clock_t my_time = clock();
  double d = find_density(fmu, kT, 1e-14, nmax);
  //printf("find_density took %g seconds\n", (clock()-my_time)/double(CLOCKS_PER_SEC));
  return d;
}

double saturated_liquid(Functional f, double kT, double nmin, double nmax,
                        double vapor_liquid_ratio) {
  Functional n = EffectivePotentialToDensity();
  double n2 = nmax;
  double n1 = nmin;
  double ftry = 0;
  double ngtry = 0;
  double ntry = 0;
  double ptry = 0;
  double pgtry = 0;
  do {
    ntry = 0.5*(n1+n2);
    ftry = f(kT, -kT*log(ntry));
    ptry = pressure(f, kT, ntry);
    if (isnan(ftry) || isnan(ptry)) {
      n2 = ntry;
    } else if (ptry < 0) {
      // If it's got a negative pressure, it's definitely got too low
      // a density, so it should be the new lower bound.
      n1 = ntry;
    } else {
      //const double mutry = find_chemical_potential(f, kT, ntry);
      ngtry = coexisting_vapor_density(f, kT, ntry, vapor_liquid_ratio);
      pgtry = pressure(f, kT, ngtry);
      //printf("nl = %g\tng = %g,\t pl = %g\t pv = %g\n", ntry, ngtry,
      //       ptry, pgtry);
      if (pgtry > ptry) {
        n1 = ntry;
      } else {
        n2 = ntry;
      }
    }
  } while (fabs(n2-n1)/n2 > 1e-12 || isnan(ftry));
  //printf("pl = %g\tpv = %g\n", pressure(f, kT, (n2+n1)*0.5), pressure(f, kT, ngtry));
  return (n2+n1)*0.5;
}


void saturated_liquid_vapor(Functional f, double kT,
                            const double nmin, const double ncrit, const double nmax,
                            double *nl_ptr, double *nv_ptr, double *mu_ptr,
                            const double fraccuracy) {
  if (isnan(*nl_ptr) || *nl_ptr < ncrit || *nl_ptr > nmax) *nl_ptr = 0.5*(nmax + ncrit);
  if (isnan(*nv_ptr) || *nv_ptr < nmin || *nv_ptr > ncrit) *nv_ptr = 0.5*(nmin + ncrit);
  // Our starting guess for mu is chosen such that there should always
  // be two minima, one on each side of ncrit, provided ncrit is
  // between the two inflection points which themselves are between
  // the two minima.
  double nl_old, nv_old, mu = find_chemical_potential(f, kT, ncrit);

  // Here I compute a start value for "closeness" that yields a stop
  // point very close to fraccuracy, to minimize the number of
  // bisections that we need to do in the best-case scenario (which
  // should happen sometimes due to our good starting guesses on late
  // minimizations.
  double closest = 0.4*fraccuracy;
  while (closest < 0.5) closest = sqrt(closest);
  int i = 0;
  do {
    nl_old = *nl_ptr;
    nv_old = *nv_ptr;

    // First solve for the optimal liquid density
    double nlmin = ncrit, nlmax = nmax;
    // First we'll see if our last guess maybe was pretty good, by
    // checking progressively closer about it.
    double nl = *nl_ptr;
    for (double closeness = closest; closeness >= fraccuracy*0.25; closeness *= closeness) {
      if (*nl_ptr*(1+closeness) > nmax || *nl_ptr*(1-closeness) < ncrit) continue;
      nl = *nl_ptr * (1.0 + closeness);
      if (find_chemical_potential(f, kT, nl) > mu) {
        // Oh well, looks like the minimum isn't in this small interval.
        nlmin = nl;
        break;
      }
      // We know the liquid density is under our "high" guess, so now
      // let's try a "low" guess just under our last solution!
      nlmax = nl;
      nl = *nl_ptr * (1.0 - closeness);
      if (find_chemical_potential(f, kT, nl) < mu) {
        // Oh well, looks like the minimum isn't in this small interval.
        nlmax = nl;
        break;
      }
      nlmin = nl;
    }
    // Now we'll just use an ordinary bisection approach.
    while (nlmax - nlmin > 0.5*fraccuracy*nlmax) {
      nl = nlmin + 0.5*(nlmax - nlmin);
      if (find_chemical_potential(f, kT, nl) > mu) nlmin = nl;
      else nlmax = nl;
    }

    // First solve for the optimal vapor density
    double nvmin = nmin, nvmax = ncrit;
    // First we'll see if our last guess maybe was pretty good, by
    // checking progressively closer about it.
    double nv = *nv_ptr;
    for (double closeness = closest; closeness > fraccuracy*0.25; closeness *= closeness) {
      if (*nv_ptr*(1+closeness) > ncrit || *nv_ptr*(1-closeness) < nmin) continue;
      nv = *nv_ptr * (1.0 + closeness);
      if (find_chemical_potential(f, kT, nv) > mu) {
        // Oh well, looks like the minimum isn't in this small interval.
        nvmin = nv;
        break;
      }
      // We know the liquid density is under our "high" guess, so now
      // let's try a "low" guess just under our last solution!
      nvmax = nv;
      nv = *nv_ptr * (1.0 - closeness);
      if (find_chemical_potential(f, kT, nv) < mu) {
        // Oh well, looks like the minimum isn't in this small interval.
        nvmax = nv;
        break;
      }
      nvmin = nv;
    }
    // Now we'll just use an ordinary bisection approach for the vapor density.
    while (nvmax - nvmin > 0.5*fraccuracy*nvmax) {
      nv = nvmin + 0.5*(nvmax - nvmin);
      if (find_chemical_potential(f, kT, nv) > mu) nvmin = nv;
      else nvmax = nv;
    }
    // Now find the slope between the two points.
    double fl = f(kT, -kT*log(nl)), fv = f(kT, -kT*log(nv));
    mu = -(fl - fv)/(nl - nv);
    *nl_ptr = nl;
    *nv_ptr = nv;
    *mu_ptr = mu;
    //printf("nl = %20.15g\tnv = %20.15g\n", nl, nv);
    if (i++ > 15) {
      printf("PANICKING in saturated_liquid_vapor!\n");
      break;
    }
  } while (fabs(*nl_ptr - nl_old) > fraccuracy*nl_old ||
           fabs(*nv_ptr - nv_old) > fraccuracy*nv_old);
}

void saturated_liquid_properties(Functional f, LiquidProperties *prop) {
  double mu;
  saturated_liquid_vapor(f, prop->kT,
                         prop->vapor_density*0.1,
                         prop->critical_density,
                         prop->liquid_density*1.2,
                         &prop->liquid_density, &prop->vapor_density, &mu, 1e-12);
}

void other_equation_of_state(FILE *o, Functional f, double kT, double nmin, double nmax) {
  Functional n = EffectivePotentialToDensity();
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
    double e = f(kT, V);
    fprintf(o, "%g\t%g\t%g\n", n, p, e);
  }
}

void equation_of_state(FILE *o, Functional f, double kT, double nmin, double nmax) {
  const double factor = 1.04;
  for (double ngoal=nmin; ngoal<nmax; ngoal *= factor) {
    double V = -kT*log(ngoal);
    double p = pressure(f, kT, ngoal);
    double der = -f.derive(kT, V)*kT/ngoal;
    //printf("Got ngoal of %g but mu is %g\n", ngoal, mu);
    double e = f(kT, V);
    fprintf(o, "%g\t%g\t%g\t%g\n", ngoal, p, e, der);
  }
}
