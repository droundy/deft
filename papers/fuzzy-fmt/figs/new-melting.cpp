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
#include "vector3d.h"

//Number of points for Monte-Carlo
const long NUM_POINTS = 80;

// radius we need to integrate around a gaussian, in units of gw.
const double inclusion_radius = 4.0;

static inline double time() {
  return clock()/double(CLOCKS_PER_SEC);
}
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

//%%%%%%%%%%%NEW ENERGY FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

struct weight {
  double n_0;
  double n_1;
  double n_2;
  double n_3;
  vector3d nv_1;
  vector3d nv_2;
};

static inline double density_gaussian(double r, double gwidth, double norm) {
  return norm*exp(-r*r*(0.5/(gwidth*gwidth)));
}

static inline double find_alpha(double temp) {
  const double sigma=1;
  const double epsilon=1;
  return sigma*pow(2/(1+sqrt((temp*log(2))/epsilon)),1.0/6);
}

static inline double find_zeta(double temp) {
  const double epsilon=1;
  const double alpha = find_alpha(temp);
  return alpha/(6*sqrt(M_PI)*(sqrt(epsilon*log(2)/temp)+log(2)));
}

weight find_weights(vector3d r, vector3d rp, double temp) {
  vector3d rdiff=r-rp;
  double rdiff_magnitude=rdiff.norm();
  const double alpha = find_alpha(temp);
  const double zeta = find_zeta(temp);
  weight w;
  w.n_2=(1/(zeta*sqrt(M_PI)))*exp(-uipow((rdiff_magnitude-(alpha/2))/zeta,2));  //ASK - should these be in the if statement below as well?
  w.n_3=(1.0/2)*(1-erf((rdiff_magnitude-(alpha/2))/zeta));
  if (rdiff_magnitude > 0) {
    w.n_0=w.n_2/(4*M_PI*rdiff_magnitude*rdiff_magnitude);
    w.n_1=w.n_2/(4*M_PI*rdiff_magnitude);
    w.nv_1 = w.n_1*(rdiff/rdiff_magnitude);
    w.nv_2 = w.n_2*(rdiff/rdiff_magnitude);
  } else {
    w.n_0=0;
    w.n_1=0;
    w.nv_1 = vector3d(0,0,0);
    w.nv_2 = vector3d(0,0,0);
  }
  return w;
}

// radius_of_peak tells us how far we need to integrate away from a
// peak with width gwidth, if the temperature is T, when finding the
// weighted densities.  This is basically asking when the weighting
// functions (see find_weights above) are negligible, which thus
// depends on the alpha and zeta parameters above.
static inline double radius_of_peak(double gwidth, double T) {
  const double alpha = find_alpha(T);
  const double zeta = find_zeta(T);
  return 0.5*alpha + 3*zeta + inclusion_radius*gwidth;
}

weight find_weighted_den_aboutR(vector3d r, vector3d R, double dx, double temp,
                                double lattice_constant, double gwidth, double norm,
                                double reduced_density) {
  const vector3d lattice_vectors[3] = {
    vector3d(0,lattice_constant/2,lattice_constant/2),
    vector3d(lattice_constant/2,0,lattice_constant/2),
    vector3d(lattice_constant/2,lattice_constant/2,0),
  };

  const int inc_Ntot= (inclusion_radius*gwidth/dx) +1; //round up! Number of infinitesimal lengths along one of the lattice_vectors

  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};
  if ((r-R).norm() > radius_of_peak(gwidth, temp)) {
    return w_den_R;
  }

  const double df = dx/(lattice_constant/2);  //sets an infinitesimal length larger than dx along lattice vector
  const vector3d da1 = lattice_vectors[0]*df; //infinitesimal lattice vectors of length df
  const vector3d da2 = lattice_vectors[1]*df;
  const vector3d da3 = lattice_vectors[2]*df;
  const double dVp = da1.cross(da2).dot(da3); //volume of infinitesimal parallelpiped
  assert(dVp > 0);
  for (int l=-inc_Ntot; l<=inc_Ntot; l++) { //integrate only over infinitesimal lattice vectors within 2 x inlusion radius
    for (int m=-inc_Ntot; m<=inc_Ntot; m++) {
      for (int o=-inc_Ntot; o<=inc_Ntot; o++) {
        const vector3d rp_from_R = l*da1 + m*da2 + o*da3;
        const vector3d rp = R + rp_from_R;
        // only bother including points within the inclusion radius:
        if (rp_from_R.norm() < inclusion_radius*gwidth) {
          weight w = find_weights(r, rp, temp);
          double n_rp = density_gaussian((rp_from_R).norm(), gwidth, norm);  // want density a distance rp-R from center of Gaussian
          w_den_R.n_0 += w.n_0*n_rp*dVp; 
          w_den_R.n_1 += w.n_1*n_rp*dVp;
          w_den_R.n_2 += w.n_2*n_rp*dVp;
          w_den_R.n_3 += w.n_3*n_rp*dVp;

          w_den_R.nv_1 += w.nv_1*n_rp*dVp;
          w_den_R.nv_2 += w.nv_2*n_rp*dVp;
        }
      }
    }
  }

  return w_den_R;
}

weight find_weighted_den_aboutR_guasquad(vector3d r, vector3d R, double dx, double temp,
    double lattice_constant,
    double gwidth, double fv) {
  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};
  if ((r-R).norm() > radius_of_peak(gwidth, temp)) {
    return w_den_R;
  }
  for (int i=-1; i<3; i=i+2) {
    for (int j=-1; j<3; j=j+2) {
      for (int k=-1; k<3; k=k+2) {
        vector3d r_prime = R+gwidth*vector3d(i, j, k);  //note: GQ abscissa=sqrt(2)/2 and this times sqrt(2)*gw = gw
        weight w = find_weights(r, r_prime, temp);

        w_den_R.n_0 += .125*(1-fv)*w.n_0;
        w_den_R.n_1 += .125*(1-fv)*w.n_1;
        w_den_R.n_2 += .125*(1-fv)*w.n_2;
        w_den_R.n_3 += .125*(1-fv)*w.n_3;

        w_den_R.nv_1 += .125*(1-fv)*w.nv_1;
        w_den_R.nv_2 += .125*(1-fv)*w.nv_2;
      }
    }
  }
  return w_den_R;
}

weight find_weighted_den_aboutR_mc(vector3d r, vector3d R, double dx, double temp,
                                   double lattice_constant,
                                   double gwidth, double fv) {
  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};
  if ((r-R).norm() > radius_of_peak(gwidth, temp)) {
    return w_den_R;
  }

  for (long i=0; i<NUM_POINTS; i++) {
    vector3d dr = vector3d::ran(gwidth);
    vector3d r_prime = R + dr;
    vector3d r_prime2 = R - dr; // using an "antithetic variate" to cancel out first-order error
    weight w = find_weights(r, r_prime, temp);
    weight w2 = find_weights(r, r_prime2, temp);

    w_den_R.n_0 += (0.5/NUM_POINTS)*(1-fv)*(w.n_0 + w2.n_0);
    w_den_R.n_1 += (0.5/NUM_POINTS)*(1-fv)*(w.n_1 + w2.n_1);
    w_den_R.n_2 += (0.5/NUM_POINTS)*(1-fv)*(w.n_2 + w2.n_2);
    w_den_R.n_3 += (0.5/NUM_POINTS)*(1-fv)*(w.n_3 + w2.n_3);

    w_den_R.nv_1 += (0.5/NUM_POINTS)*(1-fv)*(w.nv_1 + w2.nv_1);
    w_den_R.nv_2 += (0.5/NUM_POINTS)*(1-fv)*(w.nv_2 + w2.nv_2);
  }
  return w_den_R;
}

weight find_weighted_den_variances_aboutR_mc(vector3d r, vector3d R, double dx, double temp,
    double lattice_constant,
    double gwidth, double fv) {
  weight avg_w_sqr = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};
  weight avg_w = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};
  if ((r-R).norm() > radius_of_peak(gwidth, temp)) {
    return {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};
  }

  for (long i=0; i<NUM_POINTS; i++) {
    vector3d dr = vector3d::ran(gwidth);
    vector3d r_prime = R + dr;
    vector3d r_prime2 = R - dr; // using an "antithetic variate" to cancel out first-order error
    weight w = find_weights(r, r_prime, temp);
    weight w2 = find_weights(r, r_prime2, temp);

    //Compute average of squared weight:
    avg_w_sqr.n_0 += (0.25/NUM_POINTS)*sqr(w.n_0 + w2.n_0);
    avg_w_sqr.n_1 += (0.25/NUM_POINTS)*sqr(w.n_1 + w2.n_1);
    avg_w_sqr.n_2 += (0.25/NUM_POINTS)*sqr(w.n_2 + w2.n_2);
    avg_w_sqr.n_3 += (0.25/NUM_POINTS)*sqr(w.n_3 + w2.n_3);

    //Compute square of average weight:
    avg_w.n_0 += (0.5/NUM_POINTS)*(w.n_0 + w2.n_0);
    avg_w.n_1 += (0.5/NUM_POINTS)*(w.n_1 + w2.n_1);
    avg_w.n_2 += (0.5/NUM_POINTS)*(w.n_2 + w2.n_2);
    avg_w.n_3 += (0.5/NUM_POINTS)*(w.n_3 + w2.n_3);
  }

  // Compute total w variance = average of squared weight - square of average weight
  weight n_var = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};
  n_var.n_0 = (1-fv)*(avg_w_sqr.n_0 - sqr(avg_w.n_0))/NUM_POINTS;
  n_var.n_1 = (1-fv)*(avg_w_sqr.n_1 - sqr(avg_w.n_1))/NUM_POINTS;
  n_var.n_2 = (1-fv)*(avg_w_sqr.n_2 - sqr(avg_w.n_2))/NUM_POINTS;
  n_var.n_3 = (1-fv)*(avg_w_sqr.n_3 - sqr(avg_w.n_3))/NUM_POINTS;

  return n_var;
}

//find_free_energy_variance(double ?) {
//  n_var=find weighted_den_variance_aboutR_mc(r, R, dx, temp, lattice_constant, gwidth, fv);
//  double f_var += uipow(ln(1-n_3),2)*n_var.n_0 + (n_2(1-cos?)/(1-n_3))*n_var.n_1+((n_1(1-cos?)/(1-n_3))-n_2*n_2/4*M_PI*(1-n_3)*(1-n_3))*n_var.n_2+(n_0/(1-n_3)-n_1*n_2*(1-cos?)/(1-n_3)*(1-n_3)-n_2*n_2/6*M_PI*uipow(1-n_3),3))*n_var.n_3;
//}


data find_energy_new(double temp, double reduced_density, double fv, double gwidth, char *data_dir, double dx, bool verbose=false) {
  printf("\n\n#Running find_energy_new with values: temp=%g, reduced_density=%g, fv=%g, gwidth=%g, dx=%g\n", temp, reduced_density, fv, gwidth, dx);  //debug
  //printf("\nCalculating many_cells...\n");
  double reduced_num_spheres = 1-fv; // number of spheres in one primitive cell based on input vacancy fraction fv
  double lattice_constant = find_lattice_constant(reduced_density, fv);
  const vector3d lattice_vectors[3] = {
    vector3d(lattice_constant/2,lattice_constant/2,0),
    vector3d(lattice_constant/2,0,lattice_constant/2),
    vector3d(0,lattice_constant/2,lattice_constant/2),
  };
  const int Nl = (lattice_constant/2)/dx; // number of infinitesimal lengths along one of the lattice vectors
  //Nl^3 is total number of infinitesimal parallelepipeds (of volume dV) in one primitive cell
  //The volume of one infinitesimal parallelepiped dV=2dx^3 with the current definition of dx="dx_proper"/2
  //where dV=(dx_proper^3)/4 just as V=(a^3)/4 is the volume of one parallelepiped with lattice_constant=a.
  //NOTE: dx_proper chops up lattice constant a, whereas the current definition of dx chops up a/2 to create
  //chunks that correspond in number to infinitesimal lengths, or chunks, along the lattice_vectors.
  //dx_proper=lattice_constant/number of chunks along lattice vector
  //dx=(lattice_constant/2)/number of chunks along lattice vector
  const double dV = uipow(lattice_constant/Nl,3)/4.0;

  int crystal_calc_option=2;  //set to 0 for crystal free energy with brute-force integration
                              //set to 1 for crystal free energy with Gaussian Quadrature (fastest)
                              //set to 2 for crystal free energy with Monte-Carlo (more accurate)
                              
  double N_crystal = 1;  //dummy value not used if not doing brute-force integration
  if (crystal_calc_option < 1) {  // N_crystal only needs to be calculated for brute-force integration

    //Find N_crystal (number of spheres in one crystal primitive cell) to normalize reduced density n(r) later
    double N_crystal=0;
    for (int i=0; i<Nl; i++) {  //integrate over one primitive cell
      for (int j=0; j<Nl; j++) {
        for (int k=0; k<Nl; k++) {
          vector3d r=lattice_vectors[0]*i/double(Nl)
                      + lattice_vectors[1]*j/double(Nl)
                      + lattice_vectors[2]*k/double(Nl);
  
          const int many_cells=2;  //Gaussians centered at lattice points in 5x5x5 primitive cells
                                   //Gaussians father away won't contriubute much 

          for (int t=-many_cells; t <=many_cells; t++) {
            for(int u=-many_cells; u<=many_cells; u++)  {
              for (int v=-many_cells; v<= many_cells; v++) {
                const vector3d R = t*lattice_vectors[0] + u*lattice_vectors[1] + v*lattice_vectors[2];
                N_crystal += density_gaussian((r-R).norm(), gwidth, 1)*dV;   //norm=1
              }
            }
          }
        }
      }
    }
  }  //end if for N_crystal calculation
  
  if (verbose) {
    printf("Integrated number of spheres in one crystal cell is %g but we want %g\n",
           N_crystal, reduced_num_spheres);
  }
  const double norm = reduced_num_spheres/N_crystal;  //normalization constant

  //Find inhomogeneous Fideal of one crystal primitive cell
  double Fideal=0;
  for (int i=0; i<Nl; i++) {  //integrate over one primitive cell
    for (int j=0; j<Nl; j++) {
      for (int k=0; k<Nl; k++) {
        vector3d r=lattice_vectors[0]*i/double(Nl)
                   + lattice_vectors[1]*j/double(Nl)
                   + lattice_vectors[2]*k/double(Nl);

        const int many_cells=2;  //Gaussians centered at lattice points in 5x5x5 primitive cells
        //Gaussians father away won't contriubute much
        double n = 0;
        const double kT = temp;
        for (int t=-many_cells; t <=many_cells; t++) {
          for(int u=-many_cells; u<=many_cells; u++)  {
            for (int v=-many_cells; v<= many_cells; v++) {
              const vector3d R = t*lattice_vectors[0] + u*lattice_vectors[1] + v*lattice_vectors[2];
              double deltar = (r-R).norm();
              n += (1-fv)*exp(-deltar*deltar*(0.5/(gwidth*gwidth)))/uipow(sqrt(2*M_PI)*gwidth,3); // FIXME check analytic norm
            }
          }
        }
        if (n > 1e-200) { // avoid underflow and n=0 issues
          // double dF = kT*n*(log(2.646476976618268e-6*n/(sqrt(kT)*kT)) - 1.0);
          // printf("n = %g  dF = %g\n", n, dF);
          Fideal += kT*n*(log(2.646476976618268e-6*n/(sqrt(kT)*kT)) - 1.0)*dV;
        }
      }
    }
  }  //End inhomogeneous Fideal calculation
  printf("ideal gas free energy = %g\n", Fideal);

  const double max_distance_considered = radius_of_peak(gwidth, temp);
  const int many_cells = 2*max_distance_considered/lattice_constant+1;
  printf("many_cells is %d based on %g vs %g\n",
         many_cells, max_distance_considered, lattice_constant);

  //Integrate over one primitive cell (a parallelepiped) to find free energy
  double cfree_energy_per_atom;
  double cfree_energy_per_vol;
  double hfree_energy_per_atom;
  double hfree_energy_per_vol;


  for (int density_option = 0; density_option <2; density_option++) { //loop on 0 only for homogeneous free energy,
    //loop on 1 only for crystal free energy
    //loop on 0 and 1 for both

    if (density_option > 0) {
      if (crystal_calc_option < 1) {
        printf("\nCalculating Crystal Free Energy by brute-force integration...\n");
      } else if (crystal_calc_option < 2 ) {
        printf("\nCalculating Crystal Free Energy with Gaussian Quadrature...\n");
      } else printf("\nCalculating Crystal Free Energy with Monte-Carlo...\n");

      double phi_1=0, phi_2=0, phi_3=0;
      double free_energy=0;

      for (int i=0; i<Nl; i++) {
        for (int j=0; j<Nl; j++) {
          for (int k=0; k<Nl; k++) {
            vector3d r=lattice_vectors[0]*i/double(Nl)
                       + lattice_vectors[1]*j/double(Nl)
                       + lattice_vectors[2]*k/double(Nl);

            double n_0=0, n_1=0, n_2=0, n_3=0;  //weighted densities  (fundamental measures)
            //double var_n_0=0, var_n_1=0, var_n_2=0, var_n_3=0;  //variances in the weighted densities  REMOVE?
            vector3d nv_1, nv_2;
            nv_1.x=0, nv_1.y=0, nv_1.z=0, nv_2.x=0, nv_2.y=0, nv_2.z=0;
            weight n_weight= {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};  //CHECK! Initialize here?
            //weight variance_n_weight= {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0)};  //REMOVE?

            for (int t=-many_cells; t <=many_cells; t++) {
              for(int u=-many_cells; u<=many_cells; u++)  {
                for (int v=-many_cells; v<= many_cells; v++) {

                  const vector3d R = t*lattice_vectors[0] + u*lattice_vectors[1] + v*lattice_vectors[2];
                  if ((R-r).norm() < max_distance_considered) {
                    if (crystal_calc_option > 1 ) {
                      n_weight=find_weighted_den_aboutR_mc(r, R, dx, temp,  //For Crystal Free Energy in real space with Monte-Carlo
                                                           lattice_constant, gwidth, fv);
                      //variance_n_weight=find_weighted_den_variances_aboutR_mc(r, R, dx, temp,   //REMOVE?
                      //                           lattice_constant, gwidth, fv);
                    } else if (crystal_calc_option > 0 ) {
                      n_weight=find_weighted_den_aboutR_guasquad(R, r, dx, temp,  //For Crystal Free Energy in real space with Gaussian Quadrature
                               lattice_constant, gwidth, fv);
                    } else {
                      n_weight=find_weighted_den_aboutR(R, r, dx, temp,     //For Crystal Free Energy in real space without Gaussian Quadrature
                                                        lattice_constant, gwidth, norm, reduced_density);
                    }

                    // printf("Am at distance %g vs %g  with n3 contribution %g\n",
                    //        (R-r).norm(), radius_of_peak(gwidth, temp), n_weight.n_3);
                    n_0 +=n_weight.n_0;
                    n_1 +=n_weight.n_1;
                    n_2 +=n_weight.n_2;
                    n_3 +=n_weight.n_3;

                    //if (crystal_calc_option > 1 ) {   //Calculate variances when use Monte-Carlo method  REMOVE?
                    //var_n_0 +=variance_n_weight.n_0;
                    //var_n_1 +=variance_n_weight.n_1;
                    //var_n_2 +=variance_n_weight.n_2;
                    //var_n_3 +=variance_n_weight.n_3;
                    //}

                    //printf("n_weight.n_3=%g\n", n_weight.n_3);  //debug
                    //printf("n_3=%g\n", n_3);  //debug
                    if (n_3 > 1) {
                      printf("ERROR: n_3 is greater than 1 (see %g)!\n", n_3);
                      printf("  position is %g %g %g\n", r.x, r.y, r.z);
                      printf("  from just R = %g %g %g gives %g\n", R.x, R.y, R.z, n_weight.n_3);
                      exit(1);
                    }
                    // if (n_weight.n_3 > 0.2)
                    //   printf("n3(%g,%g,%g) gains %g from %g %g %g  at distance %g  i.e. %d %d %d\n",
                    //          r.x, r.y, r.z, n_weight.n_3, R.x, R.y, R.z, R.norm(), t, u, v);
                    nv_1 +=n_weight.nv_1;
                    nv_2 +=n_weight.nv_2;
                  }
                }
              }
            }
            phi_1 = -n_0*log(1-n_3);
            //printf("n_0=%g, n_3=%g, 1-n_3=%g, phi_1=%g\n", n_0, n_3, 1-n_3, phi_1);  //debug
            phi_2 = (n_1*n_2 -(nv_1.x*nv_2.x + nv_1.y*nv_2.y + nv_1.z*nv_2.z))/(1-n_3);
            //printf("n_1*n_2=%g, nv_1.x*nv_2.x=%g, 1-n_3=%g\n",n_1*n_2, nv_1.x*nv_2.x, 1-n_3);  //debug
            phi_3 = ((n_2*n_2*n_2)-(3*n_2*(nv_2.x*nv_2.x + nv_2.y*nv_2.y + nv_2.z*nv_2.z)))/(24*M_PI*(1-n_3)*(1-n_3));
            //printf("phi_1=%g, phi_2=%g, phi_3=%g\n",phi_1, phi_2, phi_3);    //debug
            free_energy += temp*(phi_1 + phi_2 + phi_3)*dV;  //NOTE: temp is actually Boltzman constant times temperature
            if (isnan(free_energy)) {
              printf("free energy is a NaN!\n");
              printf("position is: %g %g %g\n", r.x, r.y, r.z);
              printf("n0 = %g\nn1 = %g\nn2=%g\nn3=%g\n", n_0, n_1, n_2, n_3);
              printf("phi1 = %g\nphi2 = %g\nphi3=%g\n", phi_1, phi_2, phi_3);
              data data_out;
              data_out.diff_free_energy_per_atom=0;
              data_out.cfree_energy_per_atom=0;
              data_out.hfree_energy_per_vol=0;
              data_out.cfree_energy_per_vol=0;
              return data_out;
            }
            //printf("free energy is now... %g\n", free_energy);   //debug
            //printf("      finished %.5f%% of the integral\n",
            //       100*((i)/double(Nl)
            //           +(j)/uipow(Nl, 2)
            //           +(k + 1)/uipow(Nl, 3)));
          }
          ////const double fraction_complete = ((i)/double(Nl) + (j + 1)/uipow(Nl, 2));
          ////const double t = time()/60/60;
          ////const double time_total = t/fraction_complete;
          ////printf("   finished %.3f%% of the integral (%g/%g hours left)\n",
          ////       100*fraction_complete, time_total - t, time_total);
        }
        //printf("finished %.1f%% of the integral\n",
        //       100*(i + 1)/double(Nl));
        //printf("free_energy so far=%g, phi_1=%g, phi_2=%g, phi_3=%g\n",free_energy, phi_1, phi_2, phi_3);
      }

      //There are 4 parallelepipeds in 1 cube; 1 atom/parallelepiped, 4 atoms/cube;
      //4*Vol_parallelepiped=Vol_cube=lattice_constant^3
      //free_energy and Fideal are over one parallelepiped with 1-fv atoms
      cfree_energy_per_atom=(Fideal + free_energy)/reduced_num_spheres; //Fideal is the total inhomogeneous ideal free energy for 1 primitive cell
      cfree_energy_per_vol=(Fideal + free_energy)*4.0/lattice_constant*lattice_constant*lattice_constant; //
      printf("Fideal = %g\n", Fideal);
      //  --->> PUT SAME CHANGES IN find_energy?

      printf("total crystal free_energy is %g, lattice_constant is %g\n", free_energy, lattice_constant);
    }  //end if density_option > 0
    if (density_option < 1) {
      printf("\nCalculating Homogeneous Free Energy analytically ...\n");
      HomogeneousSFMTFluid hf;   //note: homogeneousFE/atom does not depend on fv or gw
      hf.sigma() = 1;
      hf.epsilon() = 1;   //energy constant in the WCA fluid
      hf.kT() = temp;
      hf.n() = reduced_density;
      hf.mu() = 0;
      //Note: hf.energy() returns energy/volume

      hfree_energy_per_atom = hf.energy()/reduced_num_spheres;   // ASK! FIX!!  free energy per sphere or "atom"
      hfree_energy_per_vol = hf.energy();    // free energy per vol
      printf("homogeneous free_energy per vol is %g\n", hf.energy());
    }
  } //end for loop - density_option

  data data_out;
  data_out.diff_free_energy_per_atom=cfree_energy_per_atom-hfree_energy_per_atom;
  data_out.cfree_energy_per_atom=cfree_energy_per_atom;
  data_out.hfree_energy_per_vol=hfree_energy_per_vol;
  data_out.cfree_energy_per_vol=cfree_energy_per_vol;

  if (crystal_calc_option > 1) {
    printf("\n*Crystal free energy calculated with Monte-Carlo\n");
  } else if (crystal_calc_option > 0) {
    printf("\n*Crystal free energy calculated with Gaussian Quadrature\n");
  } else  printf("*Crystal free energy calculated with brute-force integration\n");

  printf("*Homogeneous free energy calculated analytically\n");

  printf("data_out is: homFEperatom=%g, cryFEperatom=%g\n", hfree_energy_per_atom, data_out.cfree_energy_per_atom);
  printf("data_out is: homFEpervol=%g, cryFEpervol=%g\n", data_out.hfree_energy_per_vol, data_out.cfree_energy_per_vol);
  printf("data_out is: diffperatom=%g\n", data_out.diff_free_energy_per_atom);

  int create_alldat_file = 1;  //set to 0 for no alldat file, set to 1 to create alldat file
  // that saves data for every combination of gw and fv values ran
  if (create_alldat_file > 0)  {
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
              temp, reduced_density, fv, gwidth, hfree_energy_per_atom,
              cfree_energy_per_atom, data_out.diff_free_energy_per_atom,
              lattice_constant, reduced_num_spheres);
      fclose(newmeltoutfile);
    } else {
      printf("Unable to open file %s!\n", alldat_filename);
    }
    delete[] alldat_filename;
  }

  return data_out;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%END NEW ENERGY FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


data find_energy(double temp, double reduced_density, double fv, double gwidth, char *data_dir, double dx, bool verbose=false) {
  double reduced_num_spheres = 4*(1-fv); // number of spheres in one cell based on input vacancy fraction fv
  double vacancy = 4*fv;                 //there are 4 spheres in one cell when there are no vacancies (fv=1)
  double lattice_constant = find_lattice_constant(reduced_density, fv);
  double dV = dx*dx*dx;  //volume element dV

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
    const double norm = (1-fv)/(sqrt(2*M_PI)*gwidth*sqrt(2*M_PI)*gwidth*sqrt(2*M_PI)*gwidth); // Normalized Gaussians correspond to 4 spheres/atoms for no vacancies
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
      // at a position vector (rrx[i],rry[i],rrz[i]) from each Gaussian
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
  data_out.cfree_energy_per_vol=f.energy()/(lattice_constant*lattice_constant*lattice_constant);

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
  delete[] alldat_filename;  //ASK
  return data_out;
}

//+++++++++++++Downhill Simplex Functions+++++++++++++++

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
  printf("Evaluating Downhill Simplex Free Energies...\n");
  for (int k=0; k<3; k++) {
    //data dhill_data=find_energy(temp, reduced_density, simplex_fe[k][0], simplex_fe[k][1], data_dir, dx, verbose);
    data dhill_data=find_energy_new(temp, reduced_density, simplex_fe[k][0], simplex_fe[k][1], data_dir, dx, verbose);
    if (isnan(dhill_data.cfree_energy_per_vol)) {   //If crystal energy is NaN
      simplex_fe[k][2]=1e100;  //ASK David!
    } else simplex_fe[k][2]=dhill_data.diff_free_energy_per_atom;
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

  if (reflected.gw == 0) {  //ASK DAVID
    reflected.gw = 0.0001;
  } else if (reflected.gw < 1) {
    reflected.gw = (-1)*reflected.gw;
  }
  if (reflected.fv < 0) {  //ASK DAVID
    reflected.fv = 0;
  } else if (reflected.fv > 1) {
    reflected.fv = .99;
//  } else if (reflected.fv == 1) {   //ASK DAVID  think this condition is already taken care of somewhere else
//    reflected.fv = .99;
  }

//reflected.fe=find_energy(temp, reduced_density, reflected.fv, reflected.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
  reflected.fe=find_energy_new(temp, reduced_density, reflected.fv, reflected.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
//reflected.fe=sqrt((reflected.fv*reflected.fv) + (reflected.gw*reflected.gw));  //TEST SIMPLEX
  return reflected;
}

point_fe extend_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  point_fe extended;
  extended.fv=(3/2.0)*(simplex_fe[0][0]+simplex_fe[1][0])-(2.0*simplex_fe[2][0]);
  extended.gw=(3/2.0)*(simplex_fe[0][1]+simplex_fe[1][1])-(2.0*simplex_fe[2][1]);

  if (extended.gw == 0) {  //ASK DAVID
    extended.gw = 0.0001;
  } else if (extended.gw < 1) {
    extended.gw = (-1)*extended.gw;
  }
  if (extended.fv < 0) {  //ASK DAVID
    extended.fv = 0;
  } else if (extended.fv > 1) {
    extended.fv = .99;
//  } else if (extended.fv == 1) {   //ASK DAVID  think this condition is already taken care of somewhere else
//    extended.fv = .99;
  }

//extended.fe=find_energy(temp, reduced_density, extended.fv, extended.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
  extended.fe=find_energy_new(temp, reduced_density, extended.fv, extended.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
//extended.fe=sqrt((extended.fv*extended.fv) + (extended.gw*extended.gw));  //TEST SIMPLEX
  return extended;
}

points_fe contract_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  points_fe contracted;

  printf("working with simplex:\n"); //debug
  display_simplex(   simplex_fe);   //debug

  contracted.out.fv=((3/4.0)*(simplex_fe[0][0]+simplex_fe[1][0]))-((1/2.0)*(simplex_fe[2][0]));
  contracted.out.gw=((3/4.0)*(simplex_fe[0][1]+simplex_fe[1][1]))-((1/2.0)*(simplex_fe[2][1]));

  if (contracted.out.gw == 0) {  //ASK DAVID
    contracted.out.gw = 0.0001;
  } else if (contracted.out.gw < 1) {
    contracted.out.gw = (-1)*contracted.out.gw;
  }
  if (contracted.out.fv < 0) {  //ASK DAVID
    contracted.out.fv = 0;
  } else if (contracted.out.fv > 1) {
    contracted.out.fv = .99;
//  } else if (contracted.out.fv == 1) {   //ASK DAVID  think this condition is already taken care of somewhere else
//    contracted.out.fv = .99;
  }

  printf("contracted.out.fv=%g, contracted.out.gw=%g\n", contracted.out.fv, contracted.out.gw);   //debug
//contracted.out.fe=find_energy(temp, reduced_density, contracted.out.fv, contracted.out.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
  contracted.out.fe=find_energy_new(temp, reduced_density, contracted.out.fv, contracted.out.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
//contracted.out.fe=sqrt((contracted.out.fv*contracted.out.fv) + (contracted.out.gw*contracted.out.gw));  //TEST SIMPLEX

  contracted.in.fv=((1/4.0)*(simplex_fe[0][0]+simplex_fe[1][0]))+((1/2.0)*(simplex_fe[2][0]));
  contracted.in.gw=((1/4.0)*(simplex_fe[0][1]+simplex_fe[1][1]))+((1/2.0)*(simplex_fe[2][1]));

  if (contracted.in.gw == 0) {  //ASK DAVID
    contracted.in.gw = 0.0001;
  } else if (contracted.in.gw < 1) {
    contracted.in.gw = (-1)*contracted.in.gw;
  }
  if (contracted.in.fv < 0) {  //ASK DAVID
    contracted.in.fv = 0;
  } else if (contracted.in.fv > 1) {
    contracted.in.fv = .99;
//  } else if (contracted.in.fv == 1) {   //ASK DAVID  think this condition is already taken care of somewhere else
//    contracted.in.fv = .99;
  }

  printf("contracted.in.fv=%g, contracted.in.gw=%g\n", contracted.in.fv, contracted.in.gw);   //debug
//contracted.in.fe=find_energy(temp, reduced_density, contracted.in.fv, contracted.in.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
  contracted.in.fe=find_energy_new(temp, reduced_density, contracted.in.fv, contracted.in.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
//contracted.in.fe=sqrt((contracted.in.fv*contracted.in.fv) + (contracted.in.gw*contracted.in.gw));  //TEST SIMPLEX

  return contracted;
}

points_fe shrink_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  points_fe shrunken;

  shrunken.out.fv=(1/2.0)*(simplex_fe[0][0] + simplex_fe[1][0]);   //using in/out so don't have to make another structure
  shrunken.out.gw=(1/2.0)*(simplex_fe[0][1] + simplex_fe[1][1]);

  if (shrunken.out.gw == 0) {  //ASK DAVID
    shrunken.out.gw = 0.0001;
  } else if (shrunken.out.gw < 1) {
    shrunken.out.gw = (-1)*shrunken.out.gw;
  }
  if (shrunken.out.fv < 0) {  //ASK DAVID
    shrunken.out.fv = 0;
  } else if (shrunken.out.fv > 1) {
    shrunken.out.fv = .99;
//  } else if (shrunken.out.fv == 1) {   //ASK DAVID  think this condition is already taken care of somewhere else
//    shrunken.out.fv = .99;
  }

//shrunken.out.fe=find_energy(temp, reduced_density, shrunken.out.fv, shrunken.out.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
  shrunken.out.fe=find_energy_new(temp, reduced_density, shrunken.out.fv, shrunken.out.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
//shrunken.out.fe=sqrt((shrunken.out.fv*shrunken.out.fv) + (shrunken.out.gw*shrunken.out.gw));  //TEST SIMPLEX

  shrunken.in.fv=(1/2.0)*(simplex_fe[0][0] + simplex_fe[2][0]);
  shrunken.in.gw=(1/2.0)*(simplex_fe[0][1] + simplex_fe[2][1]);

  if (shrunken.in.gw == 0) {  //ASK DAVID
    shrunken.in.gw = 0.0001;
  } else if (shrunken.in.gw < 1) {
    shrunken.in.gw = (-1)*shrunken.in.gw;
  }
  if (shrunken.in.fv < 0) {  //ASK DAVID
    shrunken.in.fv = 0;
  } else if (shrunken.in.fv > 1) {
    shrunken.in.fv = .99;
//  } else if (shrunken.in.fv == 1) {   //ASK DAVID  think this condition is already taken care of somewhere else
//    shrunken.in.fv = .99;
  }

//shrunken.in.fe=find_energy(temp, reduced_density, shrunken.in.fv, shrunken.in.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
  shrunken.in.fe=find_energy_new(temp, reduced_density, shrunken.in.fv, shrunken.in.gw, data_dir, dx, verbose).diff_free_energy_per_atom;
//shrunken.in.fe=sqrt((shrunken.in.fv*shrunken.in.fv) + (shrunken.in.gw*shrunken.in.gw));  //TEST SIMPLEX

  return shrunken;
}

//check_convergence_simplex() {
//}

void advance_simplex(double temp, double reduced_density, double simplex_fe[3][3], char *data_dir, double dx, bool verbose) {
  printf("working with simplex:\n"); //debug
  display_simplex(   simplex_fe);   //debug
  printf("reflect simplex...\n");  //debug
  point_fe reflected_point=reflect_simplex(temp, reduced_density, simplex_fe, data_dir, dx, verbose);
  printf("reflected_point.fv=%g, reflected_point.gw=%g, reflected_point.fe=%g\n", reflected_point.fv, reflected_point.gw, reflected_point.fe);  //debug
  if (simplex_fe[0][2]  < reflected_point.fe && reflected_point.fe < simplex_fe[1][2]) {
    simplex_fe[2][0]=reflected_point.fv;
    simplex_fe[2][1]=reflected_point.gw;
    simplex_fe[2][2]=reflected_point.fe;
    printf("**simplex reflected (A):\n");  //debug
    display_simplex(   simplex_fe);   //debug
    return;
  }
  if (reflected_point.fe < simplex_fe[0][2]) {
    point_fe extended_point=extend_simplex(temp, reduced_density, simplex_fe, data_dir, dx, verbose);
    printf("extended_point.fv=%g, extended_point.gw=%g, extended_point.fe=%g\n", extended_point.fv, extended_point.gw, extended_point.fe);  //debug
    if (extended_point.fe < reflected_point.fe) {
      simplex_fe[2][0]=extended_point.fv;
      simplex_fe[2][1]=extended_point.gw;
      simplex_fe[2][2]=extended_point.fe;
      printf("**simplex reflected (B) and extended:\n");  //debug
      display_simplex(   simplex_fe);   //debug
    } else {
      simplex_fe[2][0]=reflected_point.fv;
      simplex_fe[2][1]=reflected_point.gw;
      simplex_fe[2][2]=reflected_point.fe;
      printf("**simplex reflected (C) but not extended:\n");  //debug
      display_simplex(   simplex_fe);  //debug
    }
    return;
  }
  printf("contract simplex...\n");  //debug
  printf("working with simplex:\n"); //debug
  display_simplex(   simplex_fe);   //debug
  points_fe contracted_points=contract_simplex(temp, reduced_density, simplex_fe, data_dir, dx, verbose);
  printf("contracted_points.in.fv=%g, contracted_points.in.gw=%g, contracted_points.in.fe=%g\n", contracted_points.in.fv, contracted_points.in.gw, contracted_points.in.fe);  //debug
  printf("contracted_points.out.fv=%g, contracted_points.out.gw=%g, contracted_points.out.fe=%g\n", contracted_points.out.fv, contracted_points.out.gw, contracted_points.out.fe);  //debug
  point_fe better_contracted_point;
  if (contracted_points.out.fe < contracted_points.in.fe) {
    better_contracted_point=contracted_points.out;  //there's probably a problem with this!
    printf("better_contracted_point.fv=%g, better_contracted_point.gw=%g, better_contracted_point.fe=%g\n", better_contracted_point.fv, better_contracted_point.gw, better_contracted_point.fe); //debug
  }  else {
    better_contracted_point=contracted_points.in;
  }
  printf("better_contracted_point.fv=%g, better_contracted_point.gw=%g, better_contracted_point.fe=%g\n", better_contracted_point.fv, better_contracted_point.gw, better_contracted_point.fe);  //debug
  if (better_contracted_point.fe < simplex_fe[1][2]) {
    simplex_fe[2][0]=better_contracted_point.fv;
    simplex_fe[2][1]=better_contracted_point.gw;
    simplex_fe[2][2]=better_contracted_point.fe;
    printf("**simplex contracted\n");  //debug
    display_simplex(   simplex_fe);  //debug
    return;
  }
  printf("shrink simplex...\n");  //debug
  points_fe shrunken_points=shrink_simplex(temp, reduced_density, simplex_fe, data_dir, dx, verbose);
  printf("shrunken_points.in.fv=%g, shrunken_points.in.gw=%g, shrunken_points.in.fe=%g\n", shrunken_points.in.fv, shrunken_points.in.gw, shrunken_points.in.fe);  //debug
  printf("shrunken_points.out.fv=%g, shrunken_points.out.gw=%g, shrunken_points.out.fe=%g\n", shrunken_points.out.fv, shrunken_points.out.gw, shrunken_points.out.fe);  //debug
  simplex_fe[1][0]=shrunken_points.out.fv;
  simplex_fe[1][1]=shrunken_points.out.gw;
  simplex_fe[1][2]=shrunken_points.out.fe;

  simplex_fe[2][0]=shrunken_points.in.fv;
  simplex_fe[2][1]=shrunken_points.in.gw;
  simplex_fe[2][2]=shrunken_points.in.fe;

  printf("**simplex shrunken\n");  //debug
  display_simplex(   simplex_fe);  //debug
}
//+++++++++++++++END Downhill Simplex Fucntions++++++++++++++



int main(int argc, const char **argv) {
  random::seed(1);
  double reduced_density=1.0, gw=-1, fv=-1, temp=1.0; //reduced density is the homogeneous (flat) density accounting for sphere vacancies

  //double fv_start=0.0, fv_end=.99, fv_step=0.01, gw_start=0.01, gw_end=1.5, gw_step=0.1, gw_lend=0.5, gw_lstep=0.1; //default settings
  double fv_start=0, fv_end=.1, fv_step=0.01, gw_start=0.01, gw_end=0.5, gw_step=0.01, gw_lend=0.5, gw_lstep=0.01; //default settings

  double dx=0.01;        //default grid point spacing dx=dy=dz=0.01
  int verbose = false;
  int downhill = false;

  //Downhill Simplex starting guess
  //1st col = fv,  2nd col = gw,  3rd col= Diff in Free Energy  (Crystal FE-Homogeneous FE)
  //When sorted, the first row will contain the BEST value of Free Energy Diff (lowest value),
  //the middle row will contain the MIDDLE Free Energy Diff,
  //and the last row will contain the WORST value of Free Energy Diff(highest value).
  double simplex_fe[3][3] = {{0.8,  0.1,   0},
    {0.4,  0.2,   0},
    {0.2,  0.05,  0}
  };
  //double simplex_fe[3][3] = {{80, 20, 0},  //best when ordered   TEST SIMPLEX
  //                           {40, 30, 0},   //mid when ordered
  //                           {20, 10, 0}    //worst when ordered
  //};

  char *data_dir = new char[1024];
  sprintf(data_dir,"none");
  char *default_data_dir = new char[1024];
  //sprintf(default_data_dir, "crystallization/data");
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
    {"kT", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &temp, 0, "temperature", "DOUBLE"},
    {"n", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &reduced_density, 0, "reduced density", "DOUBLE"},
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

  //*****************************End POPT Setup***********************************

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

  // Downhill Simplex
  if (downhill) {
    printf("\nInitial Downhill Simplex (1st col=fv, 2nd col=gw, 3rd col=diff free energy): \n");
    display_simplex(simplex_fe);
    //printf("Evaluating Simplex free energies (takes a bit- 14sec x 3 )...\n");
    evaluate_simplex(temp, reduced_density, simplex_fe, data_dir, dx, bool(verbose));
    display_simplex(simplex_fe);
    printf("------STARTING Downhill Simplex LOOP -----\n");
    for (int i=0; i<50; i++) {
      printf("\nLOOP %i of 50 \n", i+1);  //for debug
      printf("sort simplex...");
      sort_simplex(simplex_fe);
      printf("simplex sorted:\n");
      display_simplex(simplex_fe);  //for debug
      printf("advance simplex...");
      advance_simplex(temp, reduced_density, simplex_fe, data_dir, dx, bool(verbose));
      //printf("simplex advanced\n");  //for debug
      //display_simplex(simplex_fe);   //for debug
      //check_convergence_simplex();
    }
    sort_simplex(simplex_fe);    //delete when have loop
    display_simplex(simplex_fe);  //delete when have loop
    exit(1);
    //printf("best=%g  ", simplex_fe[0][2]);
    //printf("mid=%g  ", simplex_fe[1][2]);
    //printf("worst=%g\n", simplex_fe[2][2]);
  } //End Downhill Simplex


//TEST NEW ENERGY FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //printf("reduced_density = %g, fv = %g\n", reduced_density, fv);

  //data e_data_new =find_energy_new(temp, reduced_density, fv, gw, data_dir, dx, bool(verbose));
  //printf("e_data_new is: homFEperatom=%g, cryFEperatom=%g, diffperatom=%g, homFEpervol=%g, cryFEpervol=%g\n", e_data_new.cfree_energy_per_atom-e_data_new.diff_free_energy_per_atom, e_data_new.cfree_energy_per_atom, e_data_new.diff_free_energy_per_atom, e_data_new.hfree_energy_per_vol, e_data_new.cfree_energy_per_vol);

  //return 0;  //for debug
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//TEST Weighted Density Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Use to compare a value of weighted density caluclated analytically to that computed by
//Gaussain Quadrature (fast) and Monte-Carlo (accurate) techniques

  int do_GQ_Test = 0;  //Turn Gaussian Quadrature Test on(1)/off(0)

  if (do_GQ_Test > 0) {
    printf("reduced_density = %g, fv = %g, gw = %g\n", reduced_density, fv, gw);
    double a = find_lattice_constant(reduced_density, fv);
    vector3d r = vector3d(0,0,.55);
    vector3d R = vector3d(0,0,0);
    weight w_R = find_weighted_den_aboutR_guasquad(r, R, dx, temp, a, gw, fv);
    weight w_MC = find_weighted_den_aboutR_mc(r, R, dx, temp, a, gw, fv);

    printf("\n\nreduced_density = %g, fv = %g, gw = %g  alpha=%g zeta=%g\n", reduced_density, fv, gw,
           find_alpha(temp), find_zeta(temp));

    weight w = find_weights(r, R, temp);
    printf("w_0=%g, w_1=%g, w_2=%g, w_3=%g\n", w.n_0, w.n_1, w.n_2, w.n_3);
    w = find_weights(r+vector3d(gw,0,0), R, temp);
    printf("next w_0=%g, w_1=%g, w_2=%g, w_3=%g\n\n", w.n_0, w.n_1, w.n_2, w.n_3);

    printf("\nWeighted Density TEST results for rx=%g ry=%g  rz=%g\n", r.x, r.y, r.z);
    printf("and Gaussian point Rx=%g Ry=%g  Rz=%g  :\n", R.x, R.y, R.z);
    printf("GC: n_0=%g, n_1=%g, n_2=%g, n_3=%g\n", w_R.n_0, w_R.n_1, w_R.n_2, w_R.n_3);
    printf("MC: n_0=%g, n_1=%g, n_2=%g, n_3=%g\n", w_MC.n_0, w_MC.n_1, w_MC.n_2, w_MC.n_3);
    printf("nv_1.x=%g, nv_1.y=%g, nv_1.z=%g\n", w_R.nv_1.x, w_R.nv_1.y, w_R.nv_1.z);
    printf("nv_2.x=%g, nv_2.y=%g, nv_2.z=%g\n\n", w_R.nv_2.x, w_R.nv_2.y, w_R.nv_2.z);

    //Anayltic Solution
    double alpha=find_alpha(temp);
    double zeta=find_zeta(temp);
    double norm=uipow(gw*sqrt(2*M_PI),3);
    printf("Analytic Solution:\n");
    //n_0
    double n_0_of_r_1stterm = (uipow(gw*sqrt(2*M_PI),3)/(4*M_PI*r.z*r.z*norm*zeta*sqrt(M_PI)))*exp(-uipow((r.z-(alpha/2))/zeta,2));   //Uses first Term of w1 Taylor expansion
    double n_0_of_r = n_0_of_r_1stterm;
    printf("n_0   analytic n_0 = %g  for r.z=%g (compare with quadrature %g or mc %g)\n",
           n_0_of_r, r.z, w_R.n_0, w_MC.n_0);
    //printf("      1st term of n_0= %g\n", n_0_of_r_1stterm);
    //n_1
    double n_1_of_r_1stterm = (uipow(gw*sqrt(2*M_PI),3)/(4*M_PI*r.z*norm*zeta*sqrt(M_PI)))*exp(-uipow((r.z-(alpha/2))/zeta,2));   //Uses first Term of w1 Taylor expansion
    double n_1_of_r = n_1_of_r_1stterm;
    printf("n_1   analytic n_1 = %g  for r.z=%g (compare with quadrature %g or mc %g)\n",
           n_1_of_r, r.z, w_R.n_1, w_MC.n_1);
    //printf("      1st term of n_1= %g\n", n_1_of_r_1stterm);
    //n_2
    double n_2_of_r_1stterm = (uipow(gw*sqrt(2*M_PI),3)/(norm*zeta*sqrt(M_PI)))*exp(-uipow((r.z-(alpha/2))/zeta,2));   //Uses first Term of w2 Taylor expansion
    double n_2_of_r_3rdterm = 6*sqrt(2)*M_PI*gw*gw*gw*gw*gw/(norm*zeta*zeta*r.z)*exp(-uipow((r.z-(alpha/2))/zeta,2))*(2*uipow((r.z-(alpha/2))/zeta,2)-1);
    double n_2_of_r = n_2_of_r_1stterm + n_2_of_r_3rdterm;   //2nd term (and all even terms) = zero
    printf("n_2   analytic n_2 = %g  for r.z=%g (compare with quadrature %g or mc %g)\n",
           n_2_of_r, r.z, w_R.n_2, w_MC.n_2);
    //printf("      1st term of n_2= %g   2nd term of n_2= 0   3rd term of n_2= %g\n", n_2_of_r_1stterm, n_2_of_r_3rdterm);
    //n_3
    double n_3_of_r_1stterm = (uipow(gw*sqrt(2*M_PI),3)*(1-erf((r.z-(alpha/2))/zeta))/(2*norm));   //Uses first Term of w1 Taylor expansion
    double n_3_of_r = n_3_of_r_1stterm;
    printf("n_3   analytic n_3 = %g  for r.z=%g (compare with quadrature %g or mc %g)\n",
           n_3_of_r, r.z, w_R.n_3, w_MC.n_3);
    //printf("      1st term of n_3= %g\n", n_3_of_r_1stterm);

    //printf("alpha = %g,  zeta=%g, temp=%g\n", alpha, zeta, temp);

    //Variances
    printf("\nMonte-Carlo Variances for NUM_POINTS=%li:\n", NUM_POINTS);
    weight variance_n=find_weighted_den_variances_aboutR_mc(r, R, dx, temp, a, gw, fv);
    printf("stdev in n0 = %12g vs %12g\n", sqrt(variance_n.n_0), n_0_of_r - w_MC.n_0);
    printf("stdev in n1 = %12g vs %12g\n", sqrt(variance_n.n_1), n_1_of_r - w_MC.n_1);
    printf("stdev in n2 = %12g vs %12g\n", sqrt(variance_n.n_2), n_2_of_r - w_MC.n_2);
    printf("stdev in n3 = %12g vs %12g (compare %g)\n", sqrt(variance_n.n_3), n_3_of_r - w_MC.n_3, w_MC.n_3);

    return 0;  //for debug
  }
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    const int num_to_compute = int(0.3/0.05*1/0.01);  //program progress feedback
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
        //data e_data =find_energy(temp, reduced_density, fv, gwidth, data_dir, dx, bool(verbose));
        data e_data =find_energy_new(temp, reduced_density, fv, gwidth, data_dir, dx, bool(verbose));
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
      //data e_data =find_energy(temp, reduced_density, fv, gwidth, data_dir, dx, bool(verbose));
      data e_data =find_energy_new(temp, reduced_density, fv, gwidth, data_dir, dx, bool(verbose));
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
    //find_energy(temp, reduced_density, fv, gw, data_dir, dx, true);    //no bestdataout needed for single run
    find_energy_new(temp, reduced_density, fv, gw, data_dir, dx, bool(verbose));
  }

  return 0;
}
