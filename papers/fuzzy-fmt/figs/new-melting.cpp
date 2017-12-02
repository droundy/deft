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

static inline void took(const char *name) {
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
  double diff_free_energy_per_atom;
  double cfree_energy_per_atom;
  double hfree_energy_per_vol;
  double cfree_energy_per_vol;
};

double find_lattice_constant(double reduced_density, double fv) {
  return pow(4*(1-fv)/reduced_density, 1.0/3);
}
//%%%%%%%%%%%NEW FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

struct vec {
  double x;
  double y;
  double z;
};

struct weight {
  double n_0;
  double n_1;
  double n_2;
  double n_3;
  vec nv_1;
  vec nv_2;
};

double find_ngaus(double rx, double ry, double rz, double fv, double gwidth, double lattice_constant) {
  double n;
  const double norm = (1-fv)/pow(sqrt(2*M_PI)*gwidth,3); // Normalized Gaussians correspond to 4 spheres/atoms for no vacancies
  // multiply 4 by (1-fv) to get the reduced number of spheres.
  {
    //R1: Gaussian centered at Rx=0,     Ry=0,    Rz=0
    double dist = sqrt(rx*rx + ry*ry+rz*rz);
    n = norm*exp(-0.5*dist*dist/gwidth/gwidth);
  }
  {
    //R2: Gaussian centered at Rx=a/2,   Ry=a/2,  Rz=0
    double dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                       (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                       rz*rz);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R3: Gaussian centered at Rx=-a/2,  Ry=a/2,  Rz=0
    dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                rz*rz);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R4: Gaussian centered at Rx=a/2,   Ry=-a/2, Rz=0
    dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                rz*rz);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R5: Gaussian centered at Rx=-a/2,  Ry=-a/2, Rz=0
    dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                rz*rz);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);
  }
  {
    //R6:  Gaussian centered at Rx=0,    Ry=a/2,  Rz=a/2
    double dist = sqrt((rz-lattice_constant/2)*(rz-lattice_constant/2) +
                       (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                       rx*rx);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R7:  Gaussian centered at Rx=0,    Ry=a/2,  Rz=-a/2
    dist = sqrt((rz+lattice_constant/2)*(rz+lattice_constant/2) +
                (ry-lattice_constant/2)*(ry-lattice_constant/2) +
                rx*rx);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R8:  Gaussian centered at Rx=0,    Ry=-a/2, Rz=a/2
    dist = sqrt((rz-lattice_constant/2)*(rz-lattice_constant/2) +
                (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                rx*rx);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R9:  Gaussian centered at Rx=0,    Ry=-a/2, Rz=-a/2
    dist = sqrt((rz+lattice_constant/2)*(rz+lattice_constant/2) +
                (ry+lattice_constant/2)*(ry+lattice_constant/2) +
                rx*rx);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);
  }
  {
    //R10: Gaussian centered at Rx=a/2,  Ry=0,    Rz=a/2
    double dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                       (rz-lattice_constant/2)*(rz-lattice_constant/2) +
                       ry*ry);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R11: Gaussian centered at Rx=-a/2, Ry=0,    Rz=a/2
    dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                (rz-lattice_constant/2)*(rz-lattice_constant/2) +
                ry*ry);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R12: Gaussian centered at Rx=a/2,  Ry=0,    Rz=-a/2
    dist = sqrt((rx-lattice_constant/2)*(rx-lattice_constant/2) +
                (rz+lattice_constant/2)*(rz+lattice_constant/2) +
                ry*ry);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);

    //R13: Gaussian centered at Rx=-a/2,  Ry=0,   Rz=-a/2
    dist = sqrt((rx+lattice_constant/2)*(rx+lattice_constant/2) +
                (rz+lattice_constant/2)*(rz+lattice_constant/2) +
                ry*ry);
    n += norm*exp(-0.5*dist*dist/gwidth/gwidth);
  }
  //printf("calculated ngaus is %g\n", n);
  return n;
}

weight find_weighted_densities(double rx, double ry, double rz, double sx, double sy, double sz, int Ntot, double dx, double temp, double fv, double gwidth, double N_crystal, double reduced_density) {
  double reduced_num_spheres = 4*(1-fv); // number of spheres in one cell based on input vacancy fraction fv
  double lattice_constant = find_lattice_constant(reduced_density, fv);
  weight w_den;
  for (int l=0; l<Ntot; l++) {   //integrate over ALL space (not Ntot) FIX!
          const double rxp=l*dx +sx;
          for (int m=0; m<Ntot; m++) {
            const double ryp=m*dx +sy;
            for (int o=0; o<Ntot; o++) {
              const double rzp=o*dx +sz;
              //printf("rxp = %g, ryp= %g, rzp= %g, mag rp=%g\n", rxp, ryp, rzp, rp);

              double rxdiff=rx-rxp;
              double rydiff=ry-ryp;
              double rzdiff=rz-rzp;

              double rdiff_magnitude=pow(pow(rxdiff,2)+pow(rydiff,2)+pow(rzdiff,2),0.5);
              //printf("rxdiff_magnitude= %g, rydiff_magnitude= %g, rzdiff_magnitude= %g, mag rdiff_magnitude= %g\n", rxdiff_magnitude, rydiff_magnitude, rzdiff_magnitude, rdiff_magnitude);  //debug

              const double sigma=2;
              const double Rad=sigma/pow(2, 5.0/6);
              const double alpha = sigma*pow(2/(1+pow(temp*log(2),0.5)),1.0/6);
              const double zeta = alpha/(6*pow(M_PI,0.5)*pow((log(2)/temp),0.5)+log(2));
              //printf("alpha= %g, zeta= %g\n", alpha, zeta);   //debug
              
              double w_2=(1/zeta*pow(M_PI,2))*exp(-pow(rdiff_magnitude-(alpha/2),2)/pow(zeta,2));
              double w_0=w_2/(4*pow(M_PI,2));
              double w_1=w_2/(4*M_PI*Rad);
              double w_3=(1.0/2)*(1-erf(rdiff_magnitude));
              vec wv_1, wv_2;
              if (rdiff_magnitude > 0) {
                wv_1.x=w_1*(rxdiff/rdiff_magnitude);
                wv_1.y=w_1*(rydiff/rdiff_magnitude);
                wv_1.z=w_1*(rzdiff/rdiff_magnitude);
                wv_2.x=w_2*(rxdiff/rdiff_magnitude);
                wv_2.y=w_2*(rydiff/rdiff_magnitude);
                wv_2.z=w_2*(rzdiff/rdiff_magnitude);
              } else {
                wv_1.x=0;
                wv_1.y=0;
                wv_1.z=0;
                wv_2.x=0;
                wv_2.y=0;
                wv_2.z=0;
              }

              //printf("rxdiff_magnitude=%g, rdiff_magnitude=%g\n", rxdiff_magnitude, rdiff_magnitude);   //debug
              //printf("wv_1.x=%g, wv_2.x=%g\n", wv_1.x, wv_2.x);  //debug

              const double n_den= find_ngaus(rxp, ryp, rzp, fv, gwidth, lattice_constant)
                *(reduced_num_spheres/N_crystal);
              //printf("Crystal n_den=%g\n", n_den);
              //For homogeneous, set n_den to constant? ASK!

              double dVp=pow(dx,3);  //CHANGE THIS? - ASK!

              w_den.n_0 += n_den*w_0*dVp;
              w_den.n_1 += n_den*w_1*dVp;
              w_den.n_2 += n_den*w_2*dVp;
              w_den.n_3 += n_den*w_3*dVp;

              w_den.nv_1.x += n_den*wv_1.x*dVp;
              w_den.nv_1.y += n_den*wv_1.y*dVp;
              w_den.nv_1.z += n_den*wv_1.z*dVp;
              w_den.nv_2.x += n_den*wv_2.x*dVp;
              w_den.nv_2.y += n_den*wv_2.y*dVp;
              w_den.nv_2.z += n_den*wv_2.z*dVp;

              //printf("nv_1.x=%g, nv_2.x=%g\n", nv_1.x, nv_2.x);  //debug

              //printf("n_0= %g, n_1=%g, n_2=%g, n_3=%g\n", n_0, n_1, n_2, n_3);  //debug
            }
          }
        }
return w_den;  
}


data find_energy_new(double temp, double reduced_density, double fv, double gwidth, char *data_dir, double dx, bool verbose=false) {
  printf("\nNew find_energy function with values: temp=%g, reduced_density=%g, fv=%g, gwidth=%g, dx=%g\n", temp, reduced_density, fv, gwidth, dx);  //debug
  double reduced_num_spheres = 4*(1-fv); // number of spheres in one cell based on input vacancy fraction fv
  double lattice_constant = find_lattice_constant(reduced_density, fv);
  
  const double dV = pow(dx,3);    //ASK!
  const int Ntot=pow(((lattice_constant/dx)+1.0),3);  //number of position vectors over one cell
  printf("Ntot is %i\n", Ntot);   //debug
  
  //Normalize n(r)
  double N_crystal=0;
  for (int i=0; i<Ntot; i++) {     //integrate over one cell
    const double rx=i*dx;
    for (int j=0; j<Ntot; j++) {
      const double ry=j*dx;
      for (int k=0; k<Ntot; k++) {
        const double rz=k*dx;
        double n_den=find_ngaus(rx, ry, rz, fv, gwidth, lattice_constant);
        N_crystal += n_den*dV;
      }
    }
  }

  if (verbose) {
    printf("Integrated number of spheres in one crystal cell is %g but we want %g\n",
           N_crystal, reduced_num_spheres);
  }

  double phi_1=0, phi_2=0, phi_3=0;
  double free_energy=0;

  for (int i=0; i<Ntot; i++) {    //integrate over one cell
    const double rx=i*dx;
    for (int j=0; j<Ntot; j++) {
      const double ry=j*dx;
      for (int k=0; k<Ntot; k++) {
        const double rz=k*dx;
        //printf("rx = %g, ry= %g, rz= %g, mag r=%g\n", rx, ry, rz, r);    //debug

        double n_0=0, n_1=0, n_2=0, n_3=0;  //weighted densities  (fundamental measures)
        vec nv_1, nv_2;
        nv_1.x=0, nv_1.y=0, nv_1.z=0, nv_2.x=0, nv_2.y=0, nv_2.z=0;

        int cell_space=3;  //set cell_space variable to an odd, positive integer greater than 1
        int num_cells=pow(cell_space,3);  //total number of cells to integrate over (for our own info - not used anywhere)
        int num_cell_shifts=cell_space/2;  //used to form vectors in shifted unit cell 
        printf("ncell_space is %i, num_cells is %i, num_cell_shifts is %i\n", cell_space, num_cells, num_cell_shifts); //debug - delete later!
        
        for (int t=-1*num_cell_shifts; t < num_cell_shifts+1; t++) {
          for (int u=-1*num_cell_shifts; u < num_cell_shifts+1; u++) {
            for (int v=-1*num_cell_shifts; v < num_cell_shifts+1; v++) {
              double sx=t*lattice_constant, sy=u*lattice_constant, sz=v*lattice_constant;
              weight n_weight=find_weighted_densities(rx, ry, rz, sx, sy, sz, Ntot, dx, temp, fv, gwidth, N_crystal, reduced_density);
              n_0 +=n_weight.n_0;
              n_1 +=n_weight.n_1;
              n_2 +=n_weight.n_2;
              n_3 +=n_weight.n_3;
              nv_1.x +=n_weight.nv_1.x;
              nv_1.y +=n_weight.nv_1.y;
              nv_1.z +=n_weight.nv_1.z;
              nv_2.x +=n_weight.nv_2.x;
              nv_2.y +=n_weight.nv_2.y;
              nv_2.z +=n_weight.nv_2.z;
            }
          }
        }
        phi_1 = -n_0*log(1-n_3);
        phi_2 = (n_1*n_2 -(nv_1.x*nv_2.x + nv_1.y*nv_2.y + nv_1.z*nv_2.z))/(1-n_3);
        //printf("n_1*n_2=%g, nv_1.x*nv_2.x=%g, 1-n_3=%g\n",n_1*n_2, nv_1.x*nv_2.x, 1-n_3);  //debug
        phi_3 = (pow(n_2,3)-(3*n_2*(nv_2.x*nv_2.x + nv_2.y*nv_2.y + nv_2.z*nv_2.z)))/24*M_PI*pow((1-n_3),2);
        //printf("phi_1=%g, phi_2=%g, phi_3=%g\n",phi_1, phi_2, phi_3);    //debug
        const double epsilon=1;
        free_energy += temp*epsilon*(phi_1 + phi_2 + phi_3)*dV;
        //printf("free energy is now... %g\n", free_energy);   //debug
      }
    }
  }

  printf("free_energy is %g\n", free_energy);

  data data_out;
  data_out.diff_free_energy_per_atom=2;
  data_out.cfree_energy_per_atom=free_energy/reduced_num_spheres;   //ASK!
  data_out.hfree_energy_per_vol=2;
  data_out.cfree_energy_per_vol=free_energy/pow(lattice_constant,3);   //ASK!

  return data_out;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%END NEW FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data find_energy(double temp, double reduced_density, double fv, double gwidth, char *data_dir, double dx, bool verbose=false) {
  double reduced_num_spheres = 4*(1-fv); // number of spheres in one cell based on input vacancy fraction fv
  double vacancy = 4*fv;                 //there are 4 spheres in one cell when there are no vacancies (fv=1)
  double lattice_constant = find_lattice_constant(reduced_density, fv);
  double dV = pow(dx,3);  //volume element dV

  HomogeneousSFMTFluid hf;
  hf.sigma() = 1;
  hf.epsilon() = 1;   //energy constant in the WCA fluid
  hf.kT() = temp;
  hf.n() = reduced_density;
  hf.mu() = 0;
  //Note: hf.energy() returns energy/volume

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

  if (fv >= 1) {
    printf("Craziness: fv==1 is a meaningless scenario!");
    data data_out;
    data_out.diff_free_energy_per_atom=sqrt(-1);
    data_out.cfree_energy_per_atom=sqrt(-1);
    data_out.hfree_energy_per_vol=hf.energy();
    data_out.cfree_energy_per_vol=sqrt(-1);
    return data_out;
  }

  SFMTFluid f(lattice_constant, lattice_constant, lattice_constant, dx);
  printf("Predicted memory use: %.1f G\n", f.Nx()*f.Ny()*f.Nz()*12*8.0/1024/1024/1024);
  f.sigma() = hf.sigma();
  f.epsilon() = hf.epsilon();
  f.kT() = hf.kT();
  f.mu() = hf.mu();
  f.Vext() = 0;
  f.n() = hf.n();
  //Note: f.energy() returns energy (not energy/volume like hf.energy()!)

  double N_crystal = 0;

  {
    // This is where we set up the inhomogeneous n(r) for a Face Centered Cubic (FCC)
    const int Ntot = f.Nx()*f.Ny()*f.Nz();  //Ntot is the total number of position vectors at which the density will be calculated
    const Vector rrx = f.get_rx();          //Nx is the total number of values for rx etc...
    const Vector rry = f.get_ry();
    const Vector rrz = f.get_rz();
    const double norm = (1-fv)/pow(sqrt(2*M_PI)*gwidth,3); // Normalized Gaussians correspond to 4 spheres/atoms for no vacancies
    // multiply 4 by (1-fv) to get the reduced number of spheres.
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

  //printf("crystal free energy is %g\n", f.energy());
  double crystal_free_energy = f.energy()/reduced_num_spheres; // free energy per atom
  data data_out;
  data_out.diff_free_energy_per_atom=crystal_free_energy - homogeneous_free_energy;
  data_out.cfree_energy_per_atom=crystal_free_energy;
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
  char *alldat_filename = new char[1024];
  char *alldat_filedescriptor = new char[1024];
  sprintf(alldat_filedescriptor, "kT%5.3f_n%05.3f_fv%04.2f_gw%04.3f",
          temp, reduced_density, fv, gwidth);
  sprintf(alldat_filename, "%s/%s-alldat.dat", data_dir, alldat_filedescriptor);
  //sprintf(alldat_filename, "%s/kT%5.3f_n%05.3f_fv%04.2f_gw%04.3f-alldat.dat",
  //        data_dir, temp, reduced_density, fv, gwidth);
  printf("Create data file: %s\n", alldat_filename);

  //Create dataout file
  FILE *newmeltoutfile = fopen(alldat_filename, "w");
  if (newmeltoutfile) {
    fprintf(newmeltoutfile, "# git  version: %s\n", version_identifier());
    fprintf(newmeltoutfile, "#T\tn\tfv\tgwidth\thFreeEnergy/atom\tcFreeEnergy/atom\tFEdiff/atom\tlattice_constant\tNsph\n");
    fprintf(newmeltoutfile, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
            temp, reduced_density, fv, gwidth, homogeneous_free_energy,
            crystal_free_energy, crystal_free_energy-homogeneous_free_energy,
            lattice_constant, reduced_num_spheres);
    fclose(newmeltoutfile);
  } else {
    printf("Unable to open file %s!\n", alldat_filename);
  }
  delete[] alldat_filename;
  return data_out;
}

//+++++++++++++++++++++++++++++++Downhill Simplex+++++++++++++++++++++++++++++++

struct point_fe {
  double fv;
  double gw;
  double fe;
};

struct points_fe {
  point_fe in;
  point_fe out;
};

void display_simplex(double simplex_fe[3][3]) {
  printf("\n");
  for (int k=0; k<3; k++) {
    for(int l=0; l<3; l++) {
      printf("%g\t", simplex_fe[k][l]);
    }
    printf("\n");
  }
  printf("\n");
}


void evaluate_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  for (int k=0; k<3; k++) {
    data dhill_data=find_energy(temp, reduced_density, simplex_fe[k][0], simplex_fe[k][1], data_dir, dx, verbose);
    simplex_fe[k][2]=dhill_data.diff_free_energy_per_atom;
//simplex_fe[k][2]=sqrt((simplex_fe[k][0]*simplex_fe[k][0]) + (simplex_fe[k][1]*simplex_fe[k][1]));  //TEST SIMPLEX
    printf("simplex_fe[%i][2]=%g\n", k, simplex_fe[k][2]);
  }
}

void sort_simplex(double simplex_fe[3][3]) {
  double holdfe[3];
  for (int i =0; i < 2; ++i) {    //standard sorting algorithm
    for (int j=i+1; j < 3; j++) {
      if (simplex_fe[i][2] > simplex_fe[j][2]) {
        for (int m=0; m < 3; m++) {
          holdfe[m]=simplex_fe[i][m];
          simplex_fe[i][m]=simplex_fe[j][m];
          simplex_fe[j][m]=holdfe[m];
        }
      }
    }
  }
}

point_fe reflect_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  point_fe reflected;
  reflected.fv=simplex_fe[0][0]+simplex_fe[1][0]-simplex_fe[2][0];
  reflected.gw=simplex_fe[0][1]+simplex_fe[1][1]-simplex_fe[2][1];
  data reflect=find_energy(temp, reduced_density, reflected.fv, reflected.gw, data_dir, dx, verbose);
  reflected.fe=reflect.diff_free_energy_per_atom;
//reflected.fe=sqrt((reflected.fv*reflected.fv) + (reflected.gw*reflected.gw));  //TEST SIMPLEX
  return reflected;
}

point_fe extend_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  point_fe extended;
  extended.fv=(3/2.0)*(simplex_fe[0][0]+simplex_fe[1][0])-(2.0*simplex_fe[2][0]);
  extended.gw=(3/2.0)*(simplex_fe[0][1]+simplex_fe[1][1])-(2.0*simplex_fe[2][1]);
  data extend=find_energy(temp, reduced_density, extended.fv, extended.gw, data_dir, dx, verbose);
  extended.fe=extend.diff_free_energy_per_atom;
//extended.fe=sqrt((extended.fv*extended.fv) + (extended.gw*extended.gw));  //TEST SIMPLEX
  return extended;
}

points_fe contract_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  points_fe contracted;

  printf("working with simplex:\n"); //debug
  display_simplex(   simplex_fe);   //debug

  contracted.out.fv=((3/4.0)*(simplex_fe[0][0]+simplex_fe[1][0]))-((1/2.0)*(simplex_fe[2][0]));
  contracted.out.gw=((3/4.0)*(simplex_fe[0][1]+simplex_fe[1][1]))-((1/2.0)*(simplex_fe[2][1]));
  printf("contracted.out.fv=%g, contracted.out.gw=%g\n", contracted.out.fv, contracted.out.gw);   //debug
  data contract_out=find_energy(temp, reduced_density, contracted.out.fv, contracted.out.gw, data_dir, dx, verbose);
  contracted.out.fe=contract_out.diff_free_energy_per_atom;
//contracted.out.fe=sqrt((contracted.out.fv*contracted.out.fv) + (contracted.out.gw*contracted.out.gw));  //TEST SIMPLEX

  contracted.in.fv=((1/4.0)*(simplex_fe[0][0]+simplex_fe[1][0]))+((1/2.0)*(simplex_fe[2][0]));
  contracted.in.gw=((1/4.0)*(simplex_fe[0][1]+simplex_fe[1][1]))+((1/2.0)*(simplex_fe[2][1]));
  printf("contracted.in.fv=%g, contracted.in.gw=%g\n", contracted.in.fv, contracted.in.gw);   //debug
  data contract_in=find_energy(temp, reduced_density, contracted.in.fv, contracted.in.gw, data_dir, dx, verbose);
  contracted.in.fe=contract_in.diff_free_energy_per_atom;
//contracted.in.fe=sqrt((contracted.in.fv*contracted.in.fv) + (contracted.in.gw*contracted.in.gw));  //TEST SIMPLEX

  return contracted;
}

points_fe shrink_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  points_fe shrunken;

  shrunken.out.fv=(1/2.0)*(simplex_fe[0][0] + simplex_fe[1][0]);   //using in/out so don't have to make another structure
  shrunken.out.gw=(1/2.0)*(simplex_fe[0][1] + simplex_fe[1][1]);
  shrunken.out.fe=find_energy(temp, reduced_density, shrunken.out.fv, shrunken.out.gw,
                              data_dir, dx, verbose).diff_free_energy_per_atom;
//shrunken.out.fe=sqrt((shrunken.out.fv*shrunken.out.fv) + (shrunken.out.gw*shrunken.out.gw));  //TEST SIMPLEX

  shrunken.in.fv=(1/2.0)*(simplex_fe[0][0] + simplex_fe[2][0]);
  shrunken.in.gw=(1/2.0)*(simplex_fe[0][1] + simplex_fe[2][1]);
  shrunken.in.fe=find_energy(temp, reduced_density, shrunken.in.fv, shrunken.in.gw,
                             data_dir, dx, verbose).diff_free_energy_per_atom;
//shrunken.in.fe=sqrt((shrunken.in.fv*shrunken.in.fv) + (shrunken.in.gw*shrunken.in.gw));  //TEST SIMPLEX

  return shrunken;
}

//check_convergence_simplex() {
//}

void advance_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  printf("working with simplex:\n"); //debug
  display_simplex(   simplex_fe);   //debug
  printf("   reflect simplex\n");  //debug
  point_fe reflected_point=reflect_simplex(temp, reduced_density, simplex_fe, data_dir, dx, verbose);
  printf("   reflected_point.fv=%g, reflected_point.gw=%g, reflected_point.fe=%g\n", reflected_point.fv, reflected_point.gw, reflected_point.fe);  //debug
  if (simplex_fe[0][2]  < reflected_point.fe && reflected_point.fe < simplex_fe[1][2]) {
    simplex_fe[2][0]=reflected_point.fv;
    simplex_fe[2][1]=reflected_point.gw;
    simplex_fe[2][2]=reflected_point.fe;
    printf("   **simplex reflected A\n");  //debug
    display_simplex(   simplex_fe);   //debug
    return;
  }
  if (reflected_point.fe < simplex_fe[0][2]) {
    point_fe extended_point=extend_simplex(temp, reduced_density, simplex_fe, data_dir, dx, verbose);
    printf("   extended_point.fv=%g, extended_point.gw=%g, extended_point.fe=%g\n", extended_point.fv, extended_point.gw, extended_point.fe);  //debug
    if (extended_point.fe < reflected_point.fe) {
      simplex_fe[2][0]=extended_point.fv;
      simplex_fe[2][1]=extended_point.gw;
      simplex_fe[2][2]=extended_point.fe;
      printf("   **simplex reflected B and extended\n");  //debug
      display_simplex(   simplex_fe);   //debug
    } else {
      simplex_fe[2][0]=reflected_point.fv;
      simplex_fe[2][1]=reflected_point.gw;
      simplex_fe[2][2]=reflected_point.fe;
      printf("   **simplex reflected C but not extended\n");  //debug
      display_simplex(   simplex_fe);  //debug
    }
    return;
  }
  printf("   contract simplex\n");  //debug
  printf("working with simplex:\n"); //debug
  display_simplex(   simplex_fe);   //debug
  points_fe contracted_points=contract_simplex(temp, reduced_density, simplex_fe, data_dir, dx, verbose);
  printf("   contracted_points.in.fv=%g, contracted_points.in.gw=%g, contracted_points.in.fe=%g\n", contracted_points.in.fv, contracted_points.in.gw, contracted_points.in.fe);  //debug
  printf("   contracted_points.out.fv=%g, contracted_points.out.gw=%g, contracted_points.out.fe=%g\n", contracted_points.out.fv, contracted_points.out.gw, contracted_points.out.fe);  //debug
  point_fe better_contracted_point;
  if (contracted_points.out.fe < contracted_points.in.fe) {
    better_contracted_point=contracted_points.out;  //there's probably a problem with this!
    printf("better_contracted_point.fv=%g, better_contracted_point.gw=%g, better_contracted_point.fe=%g\n", better_contracted_point.fv, better_contracted_point.gw, better_contracted_point.fe); //debug
  }  else {
    better_contracted_point=contracted_points.in;
  }
  printf("   better_contracted_point.fv=%g, better_contracted_point.gw=%g, better_contracted_point.fe=%g\n", better_contracted_point.fv, better_contracted_point.gw, better_contracted_point.fe);  //debug
  if (better_contracted_point.fe < simplex_fe[1][2]) {
    simplex_fe[2][0]=better_contracted_point.fv;
    simplex_fe[2][1]=better_contracted_point.gw;
    simplex_fe[2][2]=better_contracted_point.fe;
    printf("   **simplex contracted\n");  //debug
    display_simplex(   simplex_fe);  //debug
    return;
  }
  printf("   shrink simplex\n");  //debug
  points_fe shrunken_points=shrink_simplex(temp, reduced_density, simplex_fe, data_dir, dx, verbose);
  printf("   shrunken_points.in.fv=%g, shrunken_points.in.gw=%g, shrunken_points.in.fe=%g\n", shrunken_points.in.fv, shrunken_points.in.gw, shrunken_points.in.fe);  //debug
  printf("   shrunken_points.out.fv=%g, shrunken_points.out.gw=%g, shrunken_points.out.fe=%g\n", shrunken_points.out.fv, shrunken_points.out.gw, shrunken_points.out.fe);  //debug
  simplex_fe[1][0]=shrunken_points.out.fv;
  simplex_fe[1][1]=shrunken_points.out.gw;
  simplex_fe[1][2]=shrunken_points.out.fe;

  simplex_fe[2][0]=shrunken_points.in.fv;
  simplex_fe[2][1]=shrunken_points.in.gw;
  simplex_fe[2][2]=shrunken_points.in.fe;

  printf("   **simplex shrunken\n");  //debug
  display_simplex(   simplex_fe);  //debug
}


//+++++++++++++++++++++++++++++++END Downhill Simplex+++++++++++++++++++++++++++



int main(int argc, const char **argv) {
  double reduced_density=1.0, gw=-1, fv=-1, temp=1.0; //reduced density is the homogeneous (flat) density accounting for sphere vacancies
  double fv_start=0.0, fv_end=.99, fv_step=0.01, gw_start=0.01, gw_end=1.5, gw_step=0.1, gw_lend=0.5, gw_lstep=0.1;
  double dx=0.01;        //grid point spacing dx=dy=dz=0.01
  int verbose = false;
  int downhill = false;

//+++++++++++++++++++++++++++++++Downhill Simplex+++++++++++++++++++++++++++++++
  double simplex_fe[3][3] = {{0.8, 0.2, 0},  //best when ordered
    {0.4, 0.3, 0},   //mid when ordered
    {0.2, 0.1, 0}    //worst when ordered
  };

//  double simplex_fe[3][3] = {{80, 20, 0},  //best when ordered   TEST SIMPLEX
//                             {40, 30, 0},   //mid when ordered
//                             {20, 10, 0}    //worst when ordered
//  };


//+++++++++++++++++++++++++++++++END Downhill Simplex+++++++++++++++++++++++++++


  char *data_dir = new char[1024];
  sprintf(data_dir,"none");
  char *default_data_dir = new char[1024];
//  sprintf(default_data_dir, "crystallization/data");
  sprintf(default_data_dir, "crystallization");


  //********************Setup POPT to get inputs from command line*******************

  // ----------------------------------------------------------------------------
  // Parse input options
  // ----------------------------------------------------------------------------

  poptOption optionsTable[] = {

    {
      "verbose", '\0', POPT_ARG_NONE, &verbose, 0,
      "Print lots of good stuff!", "BOOLEAN"
    },

    /*** FLUID PARAMETERS ***/
    {"kT", '\0', POPT_ARG_DOUBLE, &temp, 0, "temperature", "DOUBLE"},
    {"n", '\0', POPT_ARG_DOUBLE, &reduced_density, 0, "reduced density", "DOUBLE"},
    {"fv", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &fv, 0, "fraction of vacancies", "DOUBLE or -1 for loop"},
    {"gw", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw, 0, "width of Gaussian", "DOUBLE or -1 for loop (without lattice ref) or -2 loop (with lattice ref)"},

    /*** LOOPING OPTIONS ***/
    {"fvstart", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &fv_start, 0, "start fv loop at", "DOUBLE"},
    {"fvend", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &fv_end, 0, "end fv loop at", "DOUBLE"},
    {"fvstep", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &fv_step, 0, "fv loop step", "DOUBLE"},

    {"gwstart", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw_start, 0, "start gwidth loop at", "DOUBLE"},
    {"gwlend", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw_lend, 0, "end gwidth loop at lattice_constant*gw_lend", "DOUBLE"},
    {"gwlstep", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw_lstep, 0, "step by lattice_constant*gw_lstep", "DOUBLE"},
    {"gwend", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw_end, 0, "end gwidth loop at", "DOUBLE"},
    {"gwstep", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw_step, 0, "gwidth loop step", "DOUBLE"},

    /*** Downhill Simplex OPTIONS ***/
    {"dh", '\0', POPT_ARG_NONE | POPT_ARGFLAG_SHOW_DEFAULT, &downhill, 0, "Do a Downhill Simplex", "BOOLEAN"},

    /*** GRID OPTIONS ***/
    {"dx", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &dx, 0, "grid spacing dx", "DOUBLE"},

    /*** PARAMETERS DETERMINING OUTPUT FILE DIRECTORY AND NAMES ***/
    {
      "d", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &data_dir, 0,
      "Directory in which to save data", "DIRNAME"
    },

    POPT_AUTOHELP
    POPT_TABLEEND
  };

  poptContext optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
  poptSetOtherOptionHelp(optCon, "[OPTION...]\nRequired arguments: temperature (kT), "
                         "reduced density (n)");

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
    if (argv[i][1] == '-') printf("\n");
    printf("%s ", argv[i]);
  }
  printf("\n");
  printf("------------------------------------------------------------------\n\n");

//*****************************End POPT Setup**************************************

  printf("git version: %s\n", version_identifier());
  printf("\nTemperature=%g, Reduced homogeneous density=%g, Fraction of vacancies=%g, Gaussian width=%g\n", temp, reduced_density, fv, gw);

  // Create directory for data files
  if (strcmp(data_dir,"none") == 0) {
    sprintf(data_dir,"%s",default_data_dir);
    printf("\nUsing default data directory: [deft/papers/fuzzy-fmt]/%s\n", data_dir);
  } else {
    printf("\nUsing given data directory: [deft/papers/fuzzy-fmt]/%s\n", data_dir);
  }
  mkdir(data_dir, 0777);

//+++++++++++++++++++++++++++++++Downhill Simplex+++++++++++++++++++++++++++++++

// simplex_fe[3][3] = {{0.05, 0.2, 0},   //best when sorted
//                     {0.2, 0.3, 0},   //mid when sorted
//                     {0.4, 0.1, 0}};  //worst when sorted

  if (downhill) {

    display_simplex(simplex_fe);
    printf("Calculating free energies (takes a bit- 14sec x 3 )...\n");
    evaluate_simplex(temp, reduced_density, simplex_fe, data_dir, dx, bool(verbose));
    display_simplex(simplex_fe);
    printf("starting downhill simplex loop...\n");

    for (int i=0; i<50; i++) {
      printf("\nLoop %i of 50 \n", i+1);  //for debug
      printf("sort simplex\n");
      sort_simplex(simplex_fe);
      printf("simplex sorted\n");
      display_simplex(simplex_fe);  //for debug
      printf("advance simplex\n");
      advance_simplex(temp, reduced_density, simplex_fe, data_dir, dx, bool(verbose));
      printf("simplex advanced\n");
      display_simplex(simplex_fe);   //for debug
//check_convergence_simplex();
    }
    sort_simplex(simplex_fe);    //delete when have loop
    display_simplex(simplex_fe);  //delete when have loop
    exit(1);

//printf("best=%g  ", simplex_fe[0][2]);
//printf("mid=%g  ", simplex_fe[1][2]);
//printf("worst=%g\n", simplex_fe[2][2]);

  }

//+++++++++++++++++++++++++++++++END Downhill Simplex+++++++++++++++++++++++++++



//TEST NEW ENERGY FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  temp=2;
//  reduced_density=1.2;
//  fv=0.8;
//  double gwidth=0.325;

//  data e_data_new =find_energy_new(temp, reduced_density, fv, gwidth, data_dir, dx, bool(verbose));
//  printf("e_data2 is: %g, %g, %g, %g\n", e_data_new.diff_free_energy_per_atom, e_data_new.cfree_energy_per_atom, e_data_new.hfree_energy_per_vol, e_data_new.cfree_energy_per_vol);

//  return 0;  //for debug
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (fv == -1) {
    printf("fv loop variables: fv start=%g, fv_end=%g, fv step=%g\n", fv_start, fv_end, fv_step);
  }

  if (gw == -1) {
    printf("gw loop variables: gwidth start=%g, gwidth end=%g, step=%g\n", gw_start, gw_end, gw_step);
  } else if (gw == -2) {
    printf("gw loop variables: gwidth start=%g, gwidth end=lattice constant*%g, step=lattice constant*%g\n", gw_start, gw_lend, gw_lstep);
  }

  if (fv == -1) {
    double best_energy_diff = 1e100;
    double best_fv, best_gwidth, best_lattice_constant, best_cfree_energy;
    double hfree_energy_pervol, cfree_energy_pervol;
    const int num_to_compute = int(0.3/0.05*1/0.01);
    int num_computed = 0;
    for (double fv=fv_start; fv<fv_end+fv_step; fv+=fv_step) {
      double lattice_constant = find_lattice_constant(reduced_density, fv);
      printf("lattice_constant is %g\n", lattice_constant);
      if (gw == -2) {
        gw_end=lattice_constant*gw_lend;
        gw_step=lattice_constant*gw_lstep;
      } else if (gw > 0) {
        // In this case we want to do a single gwidth value, so let us
        // set the start, end and step so as to achieve that.
        gw_start = gw;
        gw_end = gw;
        gw_step = 0.1;
      }
      if (gw < 0) {
        printf ("gw is %g\n", gw);
        printf("gwend=%g, gwstep=%g   \n\n", gw_end, gw_step);
      }

      for (double gwidth=gw_start; gwidth < gw_end +0.1*gw_step; gwidth+=gw_step) {
        data e_data =find_energy(temp, reduced_density, fv, gwidth, data_dir, dx, bool(verbose));
        num_computed += 1;
        if (num_computed % (num_to_compute/100) == 0) {
          //printf("We are %.0f%% done, best_energy_diff == %g\n", 100*num_computed/double(num_to_compute),
          //       best_energy_diff);
        }
        if (e_data.diff_free_energy_per_atom < best_energy_diff) {
          //printf("better free energy with fv %g gwidth %g and E %g\n",
          //       fv, gwidth, e_data.diff_free_energy_per_atom);
          best_energy_diff = e_data.diff_free_energy_per_atom;
          best_cfree_energy = e_data.cfree_energy_per_atom;
          best_fv = fv;
          best_gwidth = gwidth;
          best_lattice_constant=lattice_constant;
          hfree_energy_pervol=e_data.hfree_energy_per_vol;
          cfree_energy_pervol=e_data.cfree_energy_per_vol;
        }
      }
    }
    printf("Best: fv %g  gwidth %g  Energy Difference %g\n", best_fv, best_gwidth, best_energy_diff);

    //Create bestdataout filename (to be used if we are looping)
    char *bestdat_filename = new char[1024];
    sprintf(bestdat_filename, "%s/kT%05.3f_n%05.3f_best.dat", data_dir, temp, reduced_density);

    //Create bestdataout file
    printf("Create best data file: %s\n", bestdat_filename);
    FILE *newmeltbest = fopen(bestdat_filename, "w");
    if (newmeltbest) {
      fprintf(newmeltbest, "# git version: %s\n", version_identifier());
      fprintf(newmeltbest, "#kT\tn\tvacancy_fraction\tGaussian_width\thomogeneous_energy/atom\t\tbest_crystal_free_energy/atom\tbest_energy_difference/atom\tbest_lattice_constant\thomogeneous_energy/volume\tbest_crystal_energy/volume\n");
      fprintf(newmeltbest, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t\t%g\t%g\n", temp, reduced_density, best_fv, best_gwidth,
              best_cfree_energy-best_energy_diff, best_cfree_energy, best_energy_diff,
              best_lattice_constant, hfree_energy_pervol, cfree_energy_pervol);
      fclose(newmeltbest);
    } else {
      printf("Unable to open file %s!\n", bestdat_filename);
    }
    delete[] bestdat_filename;

  } else if (gw < 0) {
    double best_energy_diff = 1e100;
    double best_fv, best_gwidth, best_lattice_constant, best_cfree_energy;
    double hfree_energy_pervol, cfree_energy_pervol;
    double lattice_constant = find_lattice_constant(reduced_density, fv);
    printf("lattice_constant is %g\n", lattice_constant);
    if (gw == -2) {
      gw_end=lattice_constant*gw_lend;
      gw_step=lattice_constant*gw_lstep;
    }
    printf("gw is %g\n", gw);
    printf ("gwend=%g, gwstep=%g   \n\n", gw_end, gw_step);
    for (double gwidth=gw_start; gwidth < gw_end + 0.1*gw_step; gwidth+=gw_step) {
      data e_data =find_energy(temp, reduced_density, fv, gwidth, data_dir, dx, bool(verbose));
      if (e_data.diff_free_energy_per_atom < best_energy_diff) {
        best_energy_diff = e_data.diff_free_energy_per_atom;
        best_cfree_energy = e_data.cfree_energy_per_atom;
        best_fv = fv;
        best_gwidth = gwidth;
        best_lattice_constant=lattice_constant;
        hfree_energy_pervol=e_data.hfree_energy_per_vol;
        cfree_energy_pervol=e_data.cfree_energy_per_vol;
      }
    }
    printf("For fv %g, Best: gwidth %g  energy Difference %g\n", best_fv, best_gwidth, best_energy_diff);

    //Create bestdataout filename (to be used if we are looping)
    char *bestdat_filename = new char[1024];
    sprintf(bestdat_filename, "%s/kT%05.3f_n%05.3f_best.dat", data_dir, temp, reduced_density);

    //Create bestdataout file
    printf("Create best data file: %s\n", bestdat_filename);
    FILE *newmeltbest = fopen(bestdat_filename, "w");
    if (newmeltbest) {
      fprintf(newmeltbest, "# git version: %s\n", version_identifier());
      fprintf(newmeltbest, "#kT\tn\tvacancy_fraction\tGaussian_width\thomogeneous_energy/atom\t\tbest_crystal_free_energy/atom\tbest_energy_difference/atom\tbest_lattice_constant\thomogeneous_energy/volume\tbest_crystal_energy/volume\n");
      fprintf(newmeltbest, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t\t%g\t%g\n",
              temp, reduced_density, best_fv, best_gwidth, best_cfree_energy-best_energy_diff, best_cfree_energy, best_energy_diff, best_lattice_constant, hfree_energy_pervol, cfree_energy_pervol);
      fclose(newmeltbest);
    } else {
      printf("Unable to open file %s!\n", bestdat_filename);
    }
    delete[] bestdat_filename;

  } else {
    find_energy(temp, reduced_density, fv, gw, data_dir, dx, true);
  }

  return 0;
}
