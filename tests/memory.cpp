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

#include <sched.h> // linux-specific, for sched_setaffinity

#include <sys/time.h>
#include <sys/resource.h>

#include <stdio.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include "OptimizedFunctionals.h"
#include "ContactDensity.h"
#include "LineMinimizer.h"
#include "equation-of-state.h"
#include "utilities.h"

int retval = 0;
int numoops = 0;
static double majorfaults = 0.0;
double get_time() {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  majorfaults = usage.ru_majflt/1024.0/4.0;
  return usage.ru_utime.tv_sec + 1e-6*usage.ru_utime.tv_usec;
  // return clock()/double(CLOCKS_PER_SEC);
}

double last_time = get_time();

char *hn = new char[80+1];

double check_peak(const char *name, const char *name2, FILE *out,
                  double peakmin, double peakmax,
                  double cpu = 0, double cpurat = 1.6) {
  printf("===> Testing %s of %s <===\n", name, name2);
  
  double cputime = get_time() - last_time;
  double peak = peak_memory()/1024.0/1024;
  if (cpu)
    printf("CPU time is %g s (%.0f%%, with memory use %.0f M)\n", cputime, (cputime - cpu)/cpu*100, peak);
  else
    printf("CPU time is %g s (with memory use %.0f M)\n", cputime, peak);
  if (peak < peakmin) {
    printf("FAIL: Peak memory use of %s %s should be at least %g (but it's %g)!\n", name, name2, peakmin, peak);
    retval++;
  }
  if (peak > peakmax) {
    printf("FAIL: Peak memory use of %s %s should be under %g (but it's %g)!\n", name, name2, peakmax, peak);
    retval++;
  }
  if (cpu) {
    if (cputime < cpu/cpurat) {
      printf("OOPS: CPU time of %s %s should be at least %g!\n", name, name2, cpu/cpurat);
      //retval++;
      numoops++;
    }
    if (cputime > cpu*cpurat) {
      printf("OOPS: CPU time of %s %s should be under %g!\n", name, name2, cpu*cpurat);
      //retval++;
      numoops++;
    }
  }
  if (out) {
    fprintf(out, "%g\t%g\n", peak, cputime);
    fflush(out);
  }
  reset_peak_memory();
  last_time = get_time();
  return cputime;
}

const double R = 2.7;
const double rcav = R+R; // 11.8*R+R;

double notincavity(Cartesian r) {
  const double rad2 = r.dot(r);
  if (rad2 < rcav*rcav) {
    return 0;
  } else {
    return 1;
  }
}

double incavity(Cartesian r) {
  return 1 - notincavity(r);
}

void check_a_functional(const char *name, Functional f, const Grid &x) {
  const double kT = water_prop.kT; // room temperature in Hartree

  printf("\n***********");
  for (unsigned i=0;i<strlen(name) + 4;i++) printf("*");
  printf("\n* Working on %s *\n", name);
  for (unsigned i=0;i<strlen(name) + 4;i++) printf("*");
  printf("***********\n\n");
  fflush(stdout);

  double memE, cpuE, memG, cpuG, memP, cpuP, memPonly, cpuPonly;

  FILE *out;
  {
    char *fname = new char[1024];
    snprintf(fname, 1024, "tests/bench/good/%s.%s", name, hn);
    FILE *good = fopen(fname, "r");
    bool nocpu = false;
    if (majorfaults) {
      printf("Refusing to test cpu times, since there is swapping going on! (%g M)\n", majorfaults);
      nocpu = true;
    }
    if (!good || fscanf(good, " mem cpu %lg %lg %lg %lg %lg %lg %lg %lg",
                        &memE, &cpuE,
                        &memG, &cpuG,
                        &memP, &cpuP,
                        &memPonly, &cpuPonly) != 8) {
      printf("Unable to open file %s, will not check cpu times!\n", fname);
      nocpu = true;
    } else {
      fclose(good);
    }
    snprintf(fname, 1024, "tests/bench/good/%s.kipu", name);
    good = fopen(fname, "r");
    if (!good) {
      printf("Unable to open file %s, so I have to fail!\n", fname);
      exit(1);
    }
    if (fscanf(good, " mem cpu %lg %*g %lg %*g %lg %*g %lg %*g",
               &memE, &memG, &memP, &memPonly) != 4) {
      printf("Unable to read file %s\n", fname);
      exit(1);
    }
    fclose(good);
    if (nocpu) {
      cpuE = cpuG = cpuP = cpuPonly = 0;
    }
    snprintf(fname, 1024, "tests/bench/%s.%s", name, hn);
    out = fopen(fname, "w");
    if (!out) {
      printf("Unable to create file %s\n", fname);
      exit(1);
    }
    fprintf(out, "mem\tcpu\n");
    delete[] fname;
  }

  reset_peak_memory();
  last_time = get_time();

  f.integral(kT, x);
  //printf("\n\nEnergy of %s is %g\n", name, f.integral(x));
  check_peak("Energy", name, out, memE-0.1, memE+0.1, cpuE);

  Grid mygrad(x);
  mygrad.setZero();
  f.integralgrad(kT, x, &mygrad);
  //printf("Grad of %s is: %g\n", name, out, mygrad.norm());

  check_peak("Gradient", name, out, memG-0.1, memG+0.1, cpuG);
  
  {
    Grid mypgrad(x);
    mygrad.setZero();
    mypgrad.setZero();
    f.integralgrad(kT, x, &mygrad, &mypgrad);
  }
  check_peak("Gradient and preconditioned gradient", name, out, memP-0.1, memP+0.1, cpuP);

  f.integralpgrad(kT, x, &mygrad);
  check_peak("Preconditioned gradient", name, out, memPonly-0.1, memPonly+0.1, cpuPonly);

  fclose(out);
}

int main(int, char **argv) {
  {
    // Here I set this test to continue running on its current cpu.
    // This is a slightly hokey trick to try to avoid any timings
    // variation due to NUMA and the process bouncing from one CPU to
    // another.  Ideally we'd figure out a better way to decide on
    // which CPU to use for real jobs.
    cpu_set_t cpus;
    CPU_ZERO(&cpus);
    CPU_SET(sched_getcpu(), &cpus);
    int err = sched_setaffinity(0, sizeof(cpu_set_t), &cpus);
    if (err != 0) {
      printf("Error from sched_setaffinity: %d\n", err);
    }
  }
  // Let's figure out which machine we're running on (since we only
  // want to test CPU time if we have timings for this particular
  // machine).
  gethostname(hn, 80);

  const double kT = water_prop.kT; // room temperature in Hartree
  const double eta_one = 3.0/(4*M_PI*R*R*R);
  const double nliquid = 0.324*eta_one;
  Functional n = EffectivePotentialToDensity();
  const double mu = find_chemical_potential(HardSpheres(R)(n) + IdealGasOfVeff(), kT, nliquid);

  // Here we set up the lattice.
  const double rmax = rcav*2;
  Lattice lat(Cartesian(0,rmax,rmax), Cartesian(rmax,0,rmax), Cartesian(rmax,rmax,0));
  //Lattice lat(Cartesian(1.4*rmax,0,0), Cartesian(0,1.4*rmax,0), Cartesian(0,0,1.4*rmax));
  GridDescription gd(lat, 0.2);

  last_time = get_time();
  Grid external_potential(gd);
  // Do some pointless stuff so we can get some sort of gauge as to
  // how fast this CPU is, for comparison with other tests.
  for (int i=0; i<10; i++) {
    // Do this more times, to get a more consistent result...
    external_potential = external_potential.Ones(gd.NxNyNz);
    external_potential = external_potential.cwise().exp();
    external_potential = 13*external_potential + 3*external_potential.cwise().square();
    external_potential.fft(); // compute and toss the fft...
  }
  // And now let's set the external_potential up as we'd like it.
  external_potential.Set(incavity);
  external_potential *= 1e9;
  check_peak("Setting", "external potential", NULL, 7, 8);

  Grid constraint(gd);
  constraint.Set(notincavity);
  //Functional f1 = f0 + ExternalPotential(external_potential);
  Functional ff = constrain(constraint, IdealGasOfVeff() + (HardSpheres(R) + ChemicalPotential(mu))(n));
  
  Grid potential(gd, external_potential + 0.005*VectorXd::Ones(gd.NxNyNz));

  double eps = water_prop.epsilonAB;
  double kappa = water_prop.kappaAB;

  ff = OfEffectivePotential(SaftFluid2(R, eps, kappa, water_prop.epsilon_dispersion,
                                       water_prop.lambda_dispersion, water_prop.length_scaling, mu));
  check_a_functional("SaftFluid2", ff, potential);

  //ff = OfEffectivePotential(SaftFluid(R, eps, kappa, water_prop.epsilon_dispersion,
  //                                    water_prop.lambda_dispersion, water_prop.length_scaling, mu));
  //check_a_functional("SaftFluid", ff, potential);

  //ff = Association2(R, eps, kappa, water_prop.epsilon_dispersion,
  //                  water_prop.lambda_dispersion, water_prop.length_scaling);
  //check_a_functional("Association2", ff, potential);

  //ff = Association(R, eps, kappa, water_prop.epsilon_dispersion,
  //                 water_prop.lambda_dispersion, water_prop.length_scaling);
  //check_a_functional("Association", ff, potential);

  //ff = Dispersion2(R, water_prop.epsilon_dispersion,
  //                 water_prop.lambda_dispersion, water_prop.length_scaling);
  //check_a_functional("Dispersion2", ff, potential);

  //ff = Dispersion(R, water_prop.epsilon_dispersion,
  //                water_prop.lambda_dispersion, water_prop.length_scaling);
  //check_a_functional("Dispersion", ff, potential);

  ff = constrain(constraint, (HardSpheresWBnotensor(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff());
  check_a_functional("HardSpheresWBnotensor", ff, potential);

  ff = constrain(constraint, (HardSpheresNoTensor(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff());
  check_a_functional("HardSphereNoTensor", ff, potential);

  ff = constrain(constraint, (HardSpheresNoTensor2(R) + ChemicalPotential(mu))(n) + IdealGasOfVeff());
  check_a_functional("HardSpheresNoTensor2", ff, potential);

  if (numoops == 0) {
    printf("\n%s has no oopses!\n", argv[0]);
  } else {
    printf("\n%s sort of fails %d tests!\n", argv[0], numoops);
  }
  if (retval == 0) {
    printf("\n%s passes!\n", argv[0]);
  } else {
    printf("\n%s fails %d tests!\n", argv[0], retval);
    return retval;
  }
}
