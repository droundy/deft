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
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include <popt.h>
#include "new/SFMTFluidFast.h"
#include "new/HomogeneousSFMTFluidFast.h"
#include "new/Minimize.h"
#include "version-identifier.h"

static void took(const char *name) {
  static clock_t last_time = clock();
  clock_t t = clock();
  assert(name); // so it'll count as being used...
  double peak = peak_memory()/1024.0/1024;
  printf("\t\t%s took %g seconds and %g M memory\n", name, (t-last_time)/double(CLOCKS_PER_SEC), peak);
  fflush(stdout);
  last_time = t;
}

double inhomogeneity(Vector n) {
  double maxn = n[0];
  double minn = n[0];
  for (int i=0; i<n.get_size(); i++) {
    if (n[i] > maxn) maxn = n[i];
    if (n[i] < minn) minn = n[i];
  }
  return (maxn - minn)/fabs(minn);
}

struct data {
  double diff;
  double free_energy;
  double hfree_energy_per_vol;
  double cfree_energy_per_vol;
};

double find_lattice_constant(double reduced_density, double fv) {
  return pow(4*(1-fv)/reduced_density, 1.0/3);
}

data find_energy(double temp, double reduced_density, double fv, double gwidth, char *data_dir, bool verbose=false) {
  double reduced_num_spheres = 4*(1-fv); // number of spheres in one cell based on input vacancy fraction fv
  double vacancy = 4*fv;                 //there are 4 spheres in one cell when there are no vacancies (fv=1)
  double lattice_constant = find_lattice_constant(reduced_density, fv);

  HomogeneousSFMTFluid hf;
  hf.sigma() = 1;
  hf.epsilon() = 1;   //energy constant in the WCA fluid
  hf.kT() = temp;
  hf.n() = reduced_density;
  hf.mu() = 0;

  if (verbose) {
    //printf("Reduced homogeneous density= %g, fraction of vacancies= %g, Gaussian width= %g, temp= %g\n",
    //       reduced_density, fv, gwidth, temp);

    printf("Reduced number of spheres in one fluid cell is %g, vacancy is %g spheres.\n",
           reduced_num_spheres, vacancy);
    printf("lattice constant = %g\n", lattice_constant);
  }

  const double homogeneous_free_energy = hf.energy()/reduced_density; // energy per sphere
  if (verbose) {
    printf("Bulk energy per volume is %g\n", hf.energy());
    printf("Homogeneous free energy per sphere is %g\n", homogeneous_free_energy);
  }

  const double dx = 0.01;       //grid point spacing dx=dy=dz=0.01
  const double dV = pow(dx,3);  //volume element dV
  SFMTFluid f(lattice_constant, lattice_constant, lattice_constant, dx);
  f.sigma() = hf.sigma();
  f.epsilon() = hf.epsilon();
  f.kT() = hf.kT();
  f.mu() = hf.mu();
  f.Vext() = 0;
  f.n() = hf.n();

  double N_crystal = 0;

  {
    // This is where we set up the inhomogeneous n(r) for a Face Centered Cubic (FCC)
    const int Ntot = f.Nx()*f.Ny()*f.Nz();  //Ntot is the total number of position vectors at which the density will be calculated
    const Vector rrx = f.get_rx();          //Nx is the total number of values for rx etc...
    const Vector rry = f.get_ry();
    const Vector rrz = f.get_rz();
    const double norm = (1-fv)/pow(sqrt(2*M_PI)*gwidth,3); // Using normally normalized Gaussians would correspond to 4 spheres
    // so we need to multiply by (1-fv) to get the reduced number of spheres.
    Vector setn = f.n();

    for (int i=0; i<Ntot; i++) {
      const double rx = rrx[i];
      const double ry = rry[i];
      const double rz = rrz[i];
      setn[i] = 0.0*hf.n(); //sets initial density everywhere to a small value (zero)
      // The FCC cube is set up with one whole sphere in the center of the cube.
      // dist is the magnitude of vector r-vector R=square root of ((rx-Rx)^2 + (ry-Ry)^2 + (rz-Rz)^2)
      // where r is a position vector and R is a vector to the center of a sphere or Gaussian.
      // The following code calculates the contribution to the density
      // at a position vector (rrx[i],rry[i],rrz[i]) from each Guassian
      // and adds them to get the density at that position vector which
      // is then stored in setn[i].
      //NOTE! For this code to give proper results, the Gaussians must
      //have a width that is much smaller than the lattice constant so
      //that parts of the Gaussians that extend into the cube do not
      //extend out the other sides of the cube!
      {
        //R1: Gaussian centered at Rx=0,     Ry=0,    Rz=0
        double dist = sqrt(rx*rx + ry*ry+rz*rz);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);
      }
      {
        //R2: Gaussian centered at Rx=a/2,   Ry=a/2,  Rz=0
        double dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                           (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                           rz*rz);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R3: Gaussian centered at Rx=-a/2,  Ry=a/2,  Rz=0
        dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                    (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                    rz*rz);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R4: Gaussian centered at Rx=a/2,   Ry=-a/2, Rz=0
        dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                    (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                    rz*rz);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R5: Gaussian centered at Rx=-a/2,  Ry=-a/2, Rz=0
        dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                    (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                    rz*rz);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);
      }
      {
        //R6:  Gaussian centered at Rx=0,    Ry=a/2,  Rz=a/2
        double dist = sqrt((rz-lattice_constant/2)*(rz-lattice_constant/2) +
                           (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                           rx*rx);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R7:  Gaussian centered at Rx=0,    Ry=a/2,  Rz=-a/2
        dist = sqrt((rz+lattice_constant/2)*(rz+lattice_constant/2) +
                    (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                    rx*rx);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R8:  Gaussian centered at Rx=0,    Ry=-a/2, Rz=a/2
        dist = sqrt((rz-lattice_constant/2)*(rz-lattice_constant/2) +
                    (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                    rx*rx);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R9:  Gaussian centered at Rx=0,    Ry=-a/2, Rz=-a/2
        dist = sqrt((rz+lattice_constant/2)*(rz+lattice_constant/2) +
                    (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                    rx*rx);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);
      }
      {
        //R10: Gaussian centered at Rx=a/2,  Ry=0,    Rz=a/2
        double dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                           (rz-lattice_constant/2)*(rz-lattice_constant/2) +
                           ry*ry);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R11: Gaussian centered at Rx=-a/2, Ry=0,    Rz=a/2
        dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                    (rz-lattice_constant/2)*(rz-lattice_constant/2) +
                    ry*ry);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R12: Gaussian centered at Rx=a/2,  Ry=0,    Rz=-a/2
        dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                    (rz+lattice_constant/2)*(rz+lattice_constant/2) +
                    ry*ry);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);

        //R13: Gaussian centered at Rx=-a/2,  Ry=0,   Rz=-a/2
        dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                    (rz+lattice_constant/2)*(rz+lattice_constant/2) +
                    ry*ry);
        setn[i] += norm*exp(-0.5*dist*dist/gwidth/gwidth);
      }
      //Integrate n(r) computationally to check number of spheres in one cell
      N_crystal = (setn[i]*dV) + N_crystal;
    }
    if (verbose) {
      printf("Integrated number of spheres in one crystal cell is %g but we want %g\n",
             N_crystal, reduced_num_spheres);
    }
    setn = setn*(reduced_num_spheres/N_crystal);  //Normalizes setn
    if (verbose) {
      double checking_normalized_num_spheres = 0;
      for (int i=0; i<Ntot; i++) {
        checking_normalized_num_spheres += setn[i]*dV;
      }
      printf("Integrated number of spheres in one crystal cell is NOW %.16g and we want %.16g\n",
             checking_normalized_num_spheres, reduced_num_spheres);
    }
  }

//ASK David about this stuff!! ----------------
//  if (false) {
//    char *fname = new char[5000];
//    mkdir("papers/fuzzy-fmt/figs/new-data", 0777); // make sure the directory exists
//    snprintf(fname, 5000, "papers/fuzzy-fmt/figs/new-data/initial-melting-%04.2f-%04.2f-%04.2f.dat",
//             lattice_constant, reduced_density, temp);
//    FILE *o = fopen(fname, "w");
//    if (!o) {
//      fprintf(stderr, "error creating file %s\n", fname);
//      exit(1);
//    }
//    delete[] fname;
//    const int Nz = f.Nz();
//    Vector rz = f.get_rz();
//    Vector n = f.n();
//    for (int i=0; i<Nz/2; i++) {
//      fprintf(o, "%g\t%g\n", rz[i], n[i]);
//    }
//    fclose(o);
//  }
//---------------------------------------------  


  //printf("crystal free energy is %g\n", f.energy());
  double crystal_free_energy = f.energy()/reduced_num_spheres; // free energy per sphere
  data data_out;
  data_out.diff=crystal_free_energy - homogeneous_free_energy;
  data_out.free_energy=crystal_free_energy;
  data_out.hfree_energy_per_vol=hf.energy();
  data_out.cfree_energy_per_vol=f.energy()/pow(lattice_constant,3);
  if (verbose) {
    printf("Crystal free energy is %g\n", crystal_free_energy);

    f.printme("Crystal stuff!");
    if (f.energy() != f.energy()) {
      printf("FAIL!  nan for initial energy is bad!\n");
    }

    if (crystal_free_energy < homogeneous_free_energy) {
      printf("Crystal Free Energy is LOWER than the Liquid Cell Free Energy!!!\n\n");
    } else printf("TRY AGAIN!\n\n");
  }

  // Create all output data filename
  char *alldat_filename = new char[1024];   //ASK DAVID!!
  sprintf(alldat_filename, "%s/kT%g_rd%g_fv%04.2f_gw%04.3f-alldat.dat",
          data_dir, temp, reduced_density, fv, gwidth);
  printf("Create data file: %s\n", alldat_filename);
  printf("data dir =%s\n", data_dir);

  //Create dataout file
  FILE *newmeltoutfile = fopen(alldat_filename, "w");
  if (newmeltoutfile) {
    fprintf(newmeltoutfile, "# git version: %s\n", version_identifier());  
    fprintf(newmeltoutfile, "#T\tn\tfv\tgwidth\tNsph\tlat_con\tFhom\tFcry\tdiff\n");
    fprintf(newmeltoutfile, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
            temp, reduced_density, fv, gwidth, reduced_num_spheres, lattice_constant,
            homogeneous_free_energy, crystal_free_energy,
            crystal_free_energy-homogeneous_free_energy);
    fclose(newmeltoutfile);
  } else {
    printf("Unable to open file %s!\n", alldat_filename);
  }
  return data_out;
}

int main(int argc, char **argv) {
  double reduced_density, gwidth=0.3, fv=0, temp; //reduced density is the homogeneous (flat) density accounting for sphere vacancies
  
  double fv_start=0.0, fv_end=1, fv_step=0.01, gw_start=0.01, gw_end, gw_step=10;
  double dx=0.01;
  bool verbose;  // ASK?
  
  char *data_dir = new char[1024];
  sprintf(data_dir,"none");
  char *default_data_dir = new char[1024];
  sprintf(default_data_dir, "crystalization/data");   //ASK!
  char *filename = new char[1024];
  sprintf(filename, "none");
  //char *filename_suffix = new char[1024];   //ASK!
  //sprintf(filename_suffix, "none");
  
  //if (false) {  ??ASK DAVID!  WATCH OUT right now re-writes over data!!!  FIX!
    mkdir("crystalization", 0777); // make sure the directory exists
    mkdir("crystalization/data", 0777); // make sure the directory exists
    printf("made directory [deft/papers/fuzzy-fmt]crystalization/data\n");
  //}
  
  
  //********************Setup POPT to get inputs from command line*******************

  poptContext optCon;

  // ----------------------------------------------------------------------------
  // Parse input options
  // ----------------------------------------------------------------------------

  poptOption optionsTable[] = {

    /*** FLUID PARAMETERS ***/
    {"kT", '\0', POPT_ARG_DOUBLE, &temp, 0, "temperature", "DOUBLE"},
    {"rd", '\0', POPT_ARG_DOUBLE, &reduced_density, 0, "reduced density", "DOUBLE"},
    {"fv", '\0', POPT_ARG_DOUBLE, &fv, 0, "fraction of vacancies", "DOUBLE or -1 for loop"},
    {"gw", '\0', POPT_ARG_DOUBLE, &gwidth, 0, "width of Gaussian", "DOUBLE or -1 for loop"},

    /*** LOOPING OPTIONS ***/  //ASK DAVID -- if want loop with start/end/step info must also set fv=-1 and gw=-1 as well!
    {"fvstart", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &fv_start, 0, "start fv loop at", "DOUBLE"},
    {"fvend", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &fv_end, 0, "end fv loop at", "DOUBLE"},
    {"fvstep", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &fv_step, 0, "fv loop step", "DOUBLE"},
    
    {"gwstart", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw_start, 0, "start gwidth loop at", "DOUBLE"},
 //   {"gwend", '\0', POPT_ARG_DOUBLE, &gw_end, 0, "end gwidth loop at", "DOUBLE"},  //Think about what to do with this
    {"gwstep", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw_step, 0, "gwidth loop step", "DOUBLE"},
    
 //   /*** GRID OPTIONS ***/
 //   {"dx", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &dx, 0, "grid spacing dx", "DOUBLE"},  //? if include this must pass it to find_energy() !
    
    /*** PARAMETERS DETERMINING OUTPUT FILE DIRECTORY AND NAMES ***/
    {"dir", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &data_dir, 0,
    "Directory in which to save data", "DIRNAME"},
   // {"filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
   //  "Base of output file names", "STRING"},
   // {"filename-suffix", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,   //ASK!
   //  &filename_suffix, 0, "Output file name suffix", "STRING"},
    
    POPT_AUTOHELP
    POPT_TABLEEND
  };

  optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTION...]\nRequired arguments: temperature (kT), "
                         "reduced density (rd)");
                         
  int c = 0;
  // go through arguments, set them based on optionsTable
  while((c = poptGetNextOpt(optCon)) >= 0);
  if (c < -1) {
    fprintf(stderr, "\n%s: %s\n", poptBadOption(optCon, 0), poptStrerror(c));
    return 1;
  }
  poptFreeContext(optCon);    
    
  printf("------------------------------------------------------------------\n");
  printf("Running %s with parameters:\n", argv[0]);
  for (int i = 1; i < argc; i++) {
    if (argv[i][0] == '-') printf("\n");
    printf("%s ", argv[i]);
  }
  printf("\n");
  printf("------------------------------------------------------------------\n\n");

//*****************************End POPT Setup**************************************


  printf("git version: %s\n", version_identifier());
  printf("\nTemperature=%g, Reduced homogeneous density=%g, Fraction of vacancies=%g, Gaussian width=%g\n", temp, reduced_density, fv, gwidth);
   if (fv == -1) {
    printf("fv loop variables: fv start=%g, fv_end=%g, fv step=%g\n", fv_start, fv_end, fv_step);
  }
  if (gwidth == -1) {
    printf("gw loop variables: gwidth start=%g, gwidth end=lattice constant/2, step=lattice constant/%g\n", gw_start, gw_step);
  }
  
  // Set default data directory
  if (strcmp(data_dir,"none") == 0) {
    sprintf(data_dir,"%s\n",default_data_dir);  // ????ASK!
    printf("\nUsing default data directory: [deft/papers/fuzzy-fmt]/%s\n", data_dir);
  } else {
    mkdir(data_dir, 0777); 
    printf("\nUsing given data directory: [deft/papers/fuzzy-fmt]/%s\n", data_dir);  
  }

  //Create bestdataout filename (to be used if we are looping)
  char *bestdat_filename = new char[1024];
  sprintf(bestdat_filename, "%s/kT%g_rd%g_best.dat",
          data_dir, temp, reduced_density);
  
  if (fv == -1) {
    double best_energy = 1e100;
    double best_fv, best_gwidth, best_free_energy;
    double cFEpervol;
    const int num_to_compute = int(0.3/0.05*1/0.01);
    int num_computed = 0;
    //for (double fv=0; fv<1; fv+=0.01) {  //full run
    for (double fv=fv_start; fv<fv_end+fv_step; fv+=fv_step) {   //quick run
      double lattice_constant = find_lattice_constant(reduced_density, fv);
      printf("lattice_constant is %g\n", lattice_constant);
      //for (double gwidth=0.01; gwidth <= lattice_constant/2; gwidth+=lattice_constant/10) {   //full run
      for (double gwidth=0.01; gwidth <= lattice_constant/2; gwidth+=lattice_constant/gw_step) {   //full run
      //for (double gwidth=gw_start; gwidth <= gw_end+gw_step; gwidth+=gw_step) {   //quick run
        data e_data =find_energy(temp, reduced_density, fv, gwidth, data_dir);
        num_computed += 1;
        if (num_computed % (num_to_compute/100) == 0) {
          //printf("We are %.0f%% done, best_energy == %g\n", 100*num_computed/double(num_to_compute),
          //       best_energy);
        }
        if (e_data.diff < best_energy) {
          //printf("better free energy with fv %g gwidth %g and E %g\n",
          //       fv, gwidth, e_data.diff);
          best_energy = e_data.diff;
          best_free_energy = e_data.free_energy;
          best_fv = fv;
          best_gwidth = gwidth;
          cFEpervol=e_data.cfree_energy_per_vol;
        }
      }
    }
    printf("Best: fv %g  gwidth %g  Energy Difference %g\n", best_fv, best_gwidth, best_energy);

    //Create bestdataout file
    printf("Create best data file: %s\n", bestdat_filename);
    FILE *newmeltbest = fopen(bestdat_filename, "w");
    if (newmeltbest) {
      fprintf(newmeltbest, "# git version: %s\n", version_identifier());
      fprintf(newmeltbest, "#rd\tbest_crystal_energy_per_atom\tbest_energy_difference_per_atom\t\tbest_crystal_energy_per_volume\tvacancy_fraction\n");
      fprintf(newmeltbest, "%g\t%g\t%g\t%g\t%g\n",
              reduced_density, best_free_energy, best_energy, cFEpervol, best_fv, best_gwidth);
      fclose(newmeltbest);
    } else {
      printf("Unable to open file %s!\n", bestdat_filename);
    }

  } else if (gwidth == -1) {
    double best_energy = 1e100;
    double best_fv, best_gwidth, best_free_energy;
    double cFEpervol;
    double lattice_constant = find_lattice_constant(reduced_density, fv);
    printf("lattice_constant is %g\n", lattice_constant);
    //for (double gwidth=gw_start; gwidth <= gw_end+gw_step; gwidth+=gw_step) {   //quick run
    for (double gwidth=0.01; gwidth <= lattice_constant/2; gwidth+=lattice_constant/gw_step) {   //full run
      data e_data =find_energy(temp, reduced_density, fv, gwidth, data_dir);
      if (e_data.diff < best_energy) {
          best_energy = e_data.diff;
          best_free_energy = e_data.free_energy;
          best_fv = fv;
          best_gwidth = gwidth;
          cFEpervol=e_data.cfree_energy_per_vol;
        }
    }
    printf("For fv %g, Best: gwidth %g  energy Difference %g\n", best_fv, best_gwidth, best_energy);
    
    //Create bestdataout file
    printf("Create best data file: %s\n", bestdat_filename);
    FILE *newmeltbest = fopen(bestdat_filename, "w");
    if (newmeltbest) {
      fprintf(newmeltbest, "# git version: %s\n", version_identifier());
      fprintf(newmeltbest, "#rd\tbest_crystal_energy_per_atom\tbest_energy_difference_per_atom\t\tbest_crystal_energy_per_volume\tvacancy_fraction\n");
      fprintf(newmeltbest, "%g\t%g\t%g\t%g\t%g\n",
              reduced_density, best_free_energy, best_energy, cFEpervol, best_fv, best_gwidth);
      fclose(newmeltbest);
    } else {
      printf("Unable to open file %s!\n", bestdat_filename);
    }
    
  } else {
    find_energy(temp, reduced_density, fv, gwidth, data_dir, true);
  }

  return 0;
}
