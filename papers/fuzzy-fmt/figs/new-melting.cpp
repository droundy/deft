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
#include "new/SFMTFluidVeffFast.h"
#include "new/HomogeneousSFMTFluidFast.h"
#include "new/Minimize.h"
#include "version-identifier.h"
#include "vector3d.h"

#include "findxi.h"

//Number of points for Monte-Carlo
long NUM_POINTS = 800;
double MC_ERROR = 0.0;
double seed=1;
long mc_prefactor=50000;
long mc_constant=100;
int my_experiment =0;   //set to 1 to run my experiment
bool use_tensor_weight=true;

char *free_energy_output_file = 0;

// radius we need to integrate around a gaussian, in units of gw.
const double inclusion_radius = 6.0; //original MS
//const double inclusion_radius = 6.0;   //testing higher temps 2021

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
  tensor3d nm_2;
};

static inline double density_gaussian(double r, double gwidth, double norm) {
  return norm*exp(-r*r*(0.5/(gwidth*gwidth)));
}

//static inline double find_Xi(double temp) {
  //const double epsilon=1;
  //const double alpha = find_alpha(temp);
  //return alpha/(6*sqrt(M_PI)*(sqrt(epsilon*log(2)/temp)+log(2)));
//}

static inline weight find_weights_from_alpha_Xi(vector3d r, vector3d rp, double alpha, double Xi) {
  vector3d rdiff=r-rp;
  double rdiff_magnitude=rdiff.norm();
  weight w;
  //w.n_2=(1/(Xi*sqrt(M_PI)))*exp(-uipow((rdiff_magnitude - alpha/2)/Xi,2));
  w.n_2=(sqrt(2)/(Xi*sqrt(M_PI)))*exp(-uipow((rdiff_magnitude - alpha/2)/(Xi/sqrt(2)),2));
  w.n_0=w.n_2/(4*M_PI*rdiff_magnitude*rdiff_magnitude);
  w.n_1=w.n_2/(4*M_PI*rdiff_magnitude);
  w.nv_1 = w.n_1*(rdiff/rdiff_magnitude);
  w.nv_2 = w.n_2*(rdiff/rdiff_magnitude);
  w.nm_2 = w.n_2*(rdiff.outer(rdiff)/sqr(rdiff_magnitude) - identity_tensor()*(1.0/3));   // Schmidt, Santos (2012) - use!
  //w.nm_2 = w.n_2*(rdiff.outer(rdiff)/sqr(rdiff_magnitude));   // Groh & Schmidt don't use
  //w.n_3=(1.0/2)*(1-erf((rdiff_magnitude-(alpha/2))/Xi));
  w.n_3=(1.0/2)*(1-erf((rdiff_magnitude-(alpha/2))/(Xi/sqrt(2))));
  if (rdiff_magnitude == 0) {
    w.n_0=0;
    w.n_1=0;
    w.nv_1 = vector3d(0,0,0);
    w.nv_2 = vector3d(0,0,0);
    w.nm_2=zero_tensor();
  }
  return w;
}

static inline weight find_weights(vector3d r, vector3d rp, double temp) {
  return find_weights_from_alpha_Xi(r, rp, find_alpha(temp), find_Xi(temp));
}


// radius_of_peak tells us how far we need to integrate away from a
// peak with width gwidth, if the temperature is T, when finding the
// weighted densities.  This is basically asking when the weighting
// functions (see find_weights above) are negligible, which thus
// depends on the alpha and Xi parameters above.
static inline double radius_of_peak(double gwidth, double T) {
  const double alpha = find_alpha(T);
  const double Xi = find_Xi(T);
  return 0.5*alpha + 3*Xi + inclusion_radius*gwidth;
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

  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  if ((r-R).norm() > radius_of_peak(gwidth, temp)) {
    return w_den_R;
  }

  const double df = dx/(lattice_constant/2);  // fraction of lattice vector corresponding to dx
  printf("df=dx/(lattice_constant/2))=%g\n", df);
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
          w_den_R.nm_2 += w.nm_2*n_rp*dVp;
        }
      }
    }
  }

  return w_den_R;
}

weight find_weighted_den_aboutR_guasquad(vector3d r, vector3d R, double dx, double temp,  //dx is not used but keeping format
    double lattice_constant,
    double gwidth, double fv) {
  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
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
        w_den_R.nm_2 += .125*(1-fv)*w.nm_2;
      }
    }
  }
  return w_den_R;
}

weight find_weighted_den_aboutR_mc(vector3d r, vector3d R, double dx, double temp,  //dx is not used but keeping format
                                   double lattice_constant,
                                   double gwidth, double fv) {
  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  if ((r-R).norm() > radius_of_peak(gwidth, temp)) {
    return w_den_R;
  }

  for (long i=0; i<NUM_POINTS; i++) {
    vector3d dr = vector3d::ran(gwidth);  //A vector is randomly selected from a Gaussian distribution of width gwidth
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
    w_den_R.nm_2 += (0.5/NUM_POINTS)*(1-fv)*(w.nm_2 + w2.nm_2);
  }
  return w_den_R;
}

weight find_weighted_den_aboutR_mc_accurately(vector3d r, vector3d R,
    double gwidth, double fv, double alpha, double Xi) {
  weight n = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  double n3_sqr = 0;

  // On the following line, we include the ratio of gaussian peak
  // volume to the weight function volume so as to increase the odds
  // that we get a random point that overlaps with the weight
  // functions on our first try.
  long num_points = mc_constant + mc_prefactor*uipow(gwidth/(0.5*alpha), 3);   // HERE!
  // printf("Starting with num_points = %ld\n", num_points*4);
  long i=0;
  double n3_error;
  do {
    num_points *= 4;
    for (; i<num_points; i++) {
      vector3d dr = vector3d::ran(gwidth);  //A vector is randomly selected from a Gaussian distribution of width gwidth
      vector3d r_prime = R + dr;
      vector3d r_prime2 = R - dr; // using an "antithetic variate" to cancel out first-order error
      weight w = find_weights_from_alpha_Xi(r, r_prime, alpha, Xi);
      weight w2 = find_weights_from_alpha_Xi(r, r_prime2, alpha, Xi);

      n.n_0 += 0.5*(1-fv)*(w.n_0 + w2.n_0);
      n.n_1 += 0.5*(1-fv)*(w.n_1 + w2.n_1);
      n.n_2 += 0.5*(1-fv)*(w.n_2 + w2.n_2);
      const double n3_contribution = 0.5*(1-fv)*(w.n_3 + w2.n_3);
      n.n_3 += n3_contribution;
      n3_sqr += sqr(n3_contribution);

      n.nv_1 += 0.5*(1-fv)*(w.nv_1 + w2.nv_1);
      n.nv_2 += 0.5*(1-fv)*(w.nv_2 + w2.nv_2);
      n.nm_2 += 0.5*(1-fv)*(w.nm_2 + w2.nm_2);
    }
    // we only consider error in n3, because it is dimensionless and
    // pretty easy to reason about, and the others are closely
    // related.
    n3_error = sqrt(fabs(n3_sqr/num_points - sqr(n.n_3/num_points))/(num_points-1));  //Standard Error of the Mean (SEM)
  } while (n3_error > MC_ERROR);
  //} while (n3_error > MC_ERROR || (n3_error > 0.25*fabs(1-n.n_3/num_points) && n3_error < 1e-15));
  //avoids downward spiral n3_error > ... so n3_error becomes smaller which makes ... smaller and so n3_error is still bigger etc...
  n.n_0 /= num_points;
  n.n_1 /= num_points;
  n.n_2 /= num_points;
  n.n_3 /= num_points;
  n.nv_1 /= num_points;
  n.nv_2 /= num_points;
  n.nm_2 /= num_points;
  return n;
}





double report_my_error(vector3d r, vector3d R,     //Temporary- FOR DEBUG ONLY
                       double gwidth, double fv, double alpha, double Xi) {
  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  double n3_sqr = 0;

  // On the following line, we include the ratio of gaussian peak
  // volume to the weight function volume so as to increase the odds
  // that we get a random point that overlaps with the weight
  // functions on our first try.
  long num_points = mc_constant + mc_prefactor*uipow(gwidth/(0.5*alpha), 3);   // HERE!
  // printf("Starting with num_points = %ld\n", num_points*4);
  long i=0;
  double n3_error;
  do {
    num_points *= 4;
    for (; i<num_points; i++) {
      vector3d dr = vector3d::ran(gwidth);  //A vector is randomly selected from a Gaussian distribution of width gwidth
      vector3d r_prime = R + dr;
      vector3d r_prime2 = R - dr; // using an "antithetic variate" to cancel out first-order error
      weight w = find_weights_from_alpha_Xi(r, r_prime, alpha, Xi);
      weight w2 = find_weights_from_alpha_Xi(r, r_prime2, alpha, Xi);

      w_den_R.n_0 += 0.5*(1-fv)*(w.n_0 + w2.n_0);
      w_den_R.n_1 += 0.5*(1-fv)*(w.n_1 + w2.n_1);
      w_den_R.n_2 += 0.5*(1-fv)*(w.n_2 + w2.n_2);
      w_den_R.n_3 += 0.5*(1-fv)*(w.n_3 + w2.n_3);

      w_den_R.nv_1 += 0.5*(1-fv)*(w.nv_1 + w2.nv_1);
      w_den_R.nv_2 += 0.5*(1-fv)*(w.nv_2 + w2.nv_2);
      w_den_R.nm_2 += 0.5*(1-fv)*(w.nm_2 + w2.nm_2);

      n3_sqr += 0.25*(1-fv)*(1-fv)*sqr(w.n_3 + w2.n_3);
    }
    // we only consider error in n3, because it is dimensionless and
    // pretty easy to reason about, and the others are closely
    // related.
    n3_error = sqrt(fabs(n3_sqr/num_points - sqr(w_den_R.n_3/num_points))/num_points);  //Standard Error of the Mean (SEM)
  } while (n3_error > MC_ERROR || (n3_error > 0.25*fabs(1-w_den_R.n_3/num_points) && n3_error < 1e-15));  //whichever is smaller
  return n3_error;
}






long report_total_num_points(vector3d r, vector3d R,   //Temporary- FOR DEBUG ONLY
                             double gwidth, double fv, double alpha, double Xi) {
  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  double n3_sqr = 0;

  // On the following line, we include the ratio of gaussian peak
  // volume to the weight function volume so as to increase the odds
  // that we get a random point that overlaps with the weight
  // functions on our first try.
  long num_points = mc_constant + mc_prefactor*uipow(gwidth/(0.5*alpha), 3);   // HERE!
  // printf("Starting with num_points = %ld\n", num_points*4);
  long i=0;
  double n3_error;
  do {
    num_points *= 4;
    for (; i<num_points; i++) {
      vector3d dr = vector3d::ran(gwidth);  //A vector is randomly selected from a Gaussian distribution of width gwidth
      vector3d r_prime = R + dr;
      vector3d r_prime2 = R - dr; // using an "antithetic variate" to cancel out first-order error
      weight w = find_weights_from_alpha_Xi(r, r_prime, alpha, Xi);
      weight w2 = find_weights_from_alpha_Xi(r, r_prime2, alpha, Xi);

      w_den_R.n_0 += 0.5*(1-fv)*(w.n_0 + w2.n_0);
      w_den_R.n_1 += 0.5*(1-fv)*(w.n_1 + w2.n_1);
      w_den_R.n_2 += 0.5*(1-fv)*(w.n_2 + w2.n_2);
      w_den_R.n_3 += 0.5*(1-fv)*(w.n_3 + w2.n_3);

      w_den_R.nv_1 += 0.5*(1-fv)*(w.nv_1 + w2.nv_1);
      w_den_R.nv_2 += 0.5*(1-fv)*(w.nv_2 + w2.nv_2);
      w_den_R.nm_2 += 0.5*(1-fv)*(w.nm_2 + w2.nm_2);

      n3_sqr += 0.25*(1-fv)*(1-fv)*sqr(w.n_3 + w2.n_3);
    }
    // we only consider error in n3, because it is dimensionless and
    // pretty easy to reason about, and the others are closely
    // related.
    n3_error = sqrt(fabs(n3_sqr/num_points - sqr(w_den_R.n_3/num_points))/num_points);  //Standard Error of the Mean (SEM)
  } while (n3_error > MC_ERROR || (n3_error > 0.25*fabs(1-w_den_R.n_3/num_points) && n3_error < 1e-15));  //whichever is smaller
  return num_points;
}



///////////////////////////EXPERIMENT!!!///////////////////////////////
weight find_weighted_den_aboutR_mc_accurately_experiment(vector3d r, vector3d R,
    double gwidth, double fv, double alpha, double Xi) {
  weight w_den_R = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  printf("SIZE of Gaussian is 3*gwidth=%g\n", 3*gwidth);
  printf("SIZE of weight function is 2*Xi*sqrt(1/2) + alpha/2=%g\n", 2*Xi*sqrt(1/2) + alpha/2);
  double n3_sqr = 0;
  //int R0num_points;
  //vector3d rdiff = r - R;
  // printf("Starting with num_points = %ld\n", num_points*4);
  long i=0;
  double n3_error;
  int vec_count=0;
  long num_points = 5;
  //long num_points = mc_constant + mc_prefactor*uipow(gwidth/(0.5*alpha), 3);   // HERE!
  int num_zeroterms=0;
  //printf("again, num_points=%li  ", num_points);
  do {
    num_points *= 4;
    for (; i<num_points; i++) {
      //vector3d dr = vector3d::ran(gwidth);  //A vector is randomly selected from a Gaussian distribution of width gwidth
      vector3d dr;
      vector3d r_prime;
      vector3d r_prime2;
      //New experimenatal sampling technique
      //This uses a "sampling sphere" of radius |Amax| about weight functions within which points
      //generated randomly must fall in order to be included in the Monte-Carlo
      //average weight calculation. |Amax| was 3*(sqrt(Xi/2)) + alpha/2 but
      //now |Amax|=2*(sqrt(Xi/2)) + alpha/2 as the data fits this better
      vector3d A;
      if (3*gwidth > (2*Xi*sqrt(1/2) + alpha/2)) {
        do {
          dr = vector3d::ran(gwidth);  //A vector is randomly selected from a Gaussian distribution of width gwidth
          vec_count += 1;
          A=r-R-dr;   //Sample Shpere around weight function is of radius |Amax|
          printf("A.norm()=%g\n", A.norm());
        } while (A.norm() > 2*Xi*sqrt(1/2) + alpha/2);   //FIX THIS!! getting stuck in here!!! the Gaussian is mostly smaller than the weight function - turn this around!

        r_prime = R + dr;
        r_prime2 = R - dr; // using an "antithetic variate" to cancel out first-order error
      } else {
        dr = vector3d::ran(gwidth);
        r_prime = R + dr;
        r_prime2 = R - dr; // using an "antithetic variate" to cancel out first-order error
      }
      //if (R.norm() == 0) {
      //  R0num_points += 1; }
      //printf("COUNT=%i\n", vec_count);
      weight w = find_weights_from_alpha_Xi(r, r_prime, alpha, Xi);
      weight w2 = find_weights_from_alpha_Xi(r, r_prime2, alpha, Xi);
      //printf("A.norm=%g, w.n_0=%g, w.n_1=%g, w.n_2=%g, w.n_3=%g, 3*sqrt(Xi/2)=%g, Amax3=%g, Amax2=%g, Amax1=%g, alpha/2=%g, Amin1=%g, Amin2=%g \n",   //Check weighting function data:
      //        A.norm(), w.n_0, w.n_1, w.n_2, w.n_3, 3*sqrt(Xi/2), (alpha/2)+(3*sqrt(Xi/2)), (alpha/2)+(2*sqrt(Xi/2)), (alpha/2)+(1*sqrt(Xi/2)),
      //        alpha/2, (alpha/2)-(1*sqrt(Xi/2)), (alpha/2)-(2*sqrt(Xi/2)));
      w_den_R.n_0 += 0.5*(1-fv)*(w.n_0 + w2.n_0);
      w_den_R.n_1 += 0.5*(1-fv)*(w.n_1 + w2.n_1);
      w_den_R.n_2 += 0.5*(1-fv)*(w.n_2 + w2.n_2);
      w_den_R.n_3 += 0.5*(1-fv)*(w.n_3 + w2.n_3);

      w_den_R.nv_1 += 0.5*(1-fv)*(w.nv_1 + w2.nv_1);
      w_den_R.nv_2 += 0.5*(1-fv)*(w.nv_2 + w2.nv_2);
      w_den_R.nm_2 += 0.5*(1-fv)*(w.nm_2 + w2.nm_2);

      n3_sqr += 0.25*(1-fv)*(1-fv)*sqr(w.n_3 + w2.n_3);
    }
    // we only consider error in n3, because it is dimensionless and
    // pretty easy to reason about, and the others are closely
    // related.
    n3_error = sqrt(fabs(n3_sqr/num_points - sqr(w_den_R.n_3/num_points))/num_points);  //Standard Error of the Mean (SEM)
  } while (n3_error > MC_ERROR || (n3_error > 0.25*fabs(1-w_den_R.n_3/num_points) && n3_error < 1e-15));
  num_zeroterms=vec_count*(3*gwidth/(2*Xi*sqrt(1/2) + alpha/2))*(3*gwidth/(2*Xi*sqrt(1/2) + alpha/2))*(3*gwidth/(2*Xi*sqrt(1/2) + alpha/2));
  w_den_R.n_0 /= num_points + num_zeroterms;    //Added terms where weight function is zero but the Gaussian is not!
  w_den_R.n_1 /= num_points + num_zeroterms;    //Added terms where weight function is zero but the Gaussian is not!
  w_den_R.n_2 /= num_points + num_zeroterms;    //Added terms where weight function is zero but the Gaussian is not!
  w_den_R.n_3 /= num_points + num_zeroterms;    //Added terms where weight function is zero but the Gaussian is not!
  w_den_R.nv_1 /= num_points + num_zeroterms;   //Added terms where weight function is zero but the Gaussian is not!
  w_den_R.nv_2 /= num_points + num_zeroterms;   //Added terms where weight function is zero but the Gaussian is not!
  w_den_R.nm_2 /= num_points + num_zeroterms;   //Added terms where weight function is zero but the Gaussian is not!
  //vector3d rdiff=r-R;
  //printf("Total num_points is %li for |r-R|=%g\n", num_points, rdiff.norm());
  printf("vec_count=%i\n", vec_count);
  printf("3*gwidth=%g\n", 3*gwidth);
  printf("2*Xi*sqrt(1/2) + alpha/2=%g\n", 2*Xi*sqrt(1/2) + alpha/2);
  printf("3*gwidth/(2*Xi*sqrt(1/2) + alpha/2)=%g\n", 3*gwidth/(2*Xi*sqrt(1/2) + alpha/2));
  printf("num_zeroterms=%i\n", num_zeroterms);
  return w_den_R;   //the experiment takes much longer because the weight function is much larger than the gaussian up to gw=0.165 (at KT=2)
}
//////////////////////// END EXPERIMENT ///////////////////////////////



weight find_weighted_den_variances_aboutR_mc(vector3d r, vector3d R, double dx, double temp,  //dx is not used but keeping format
    double lattice_constant,
    double gwidth, double fv) {
  weight avg_w_sqr = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  weight avg_w = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  if ((r-R).norm() > radius_of_peak(gwidth, temp)) {
    return {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  }

  for (long i=0; i<NUM_POINTS; i++) {
    vector3d dr = vector3d::ran(gwidth); //A vector is randomly selected from a Gaussian distribution of width gwidth
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
  weight n_var = {0,0,0,0,vector3d(0,0,0), vector3d(0,0,0), zero_tensor()};
  n_var.n_0 = (1-fv)*(avg_w_sqr.n_0 - sqr(avg_w.n_0))/NUM_POINTS;
  n_var.n_1 = (1-fv)*(avg_w_sqr.n_1 - sqr(avg_w.n_1))/NUM_POINTS;
  n_var.n_2 = (1-fv)*(avg_w_sqr.n_2 - sqr(avg_w.n_2))/NUM_POINTS;
  n_var.n_3 = (1-fv)*(avg_w_sqr.n_3 - sqr(avg_w.n_3))/NUM_POINTS;

  return n_var;
}


data find_energy_new(double temp, double reduced_density, double fv, double gwidth, char *data_dir, double dx_input, bool verbose=false) {
  double start_time = time();
  printf("\n\n#Running find_energy_new with values: temp=%g, reduced_density=%g, fv=%g, gwidth=%g, dx=%g, mc NUM_POINTS=%li\n", temp, reduced_density, fv, gwidth, dx_input, NUM_POINTS);  //debug
  //printf("\nCalculating many_cells...\n");
  double reduced_num_spheres = 1-fv; // number of spheres in one primitive cell based on input vacancy fraction fv
  double lattice_constant = find_lattice_constant(reduced_density, fv);
  printf("lattice_constant=%g\n", lattice_constant);
  // const double cubic_cell_volume = uipow(lattice_constant, 3);




  const vector3d lattice_vectors[3] = {
    vector3d(0,lattice_constant/2,lattice_constant/2),
    vector3d(lattice_constant/2,0,lattice_constant/2),
    vector3d(lattice_constant/2,lattice_constant/2,0),
  };
  // Note: the primitive cell volume is precisely 25% of the cubic_cell_volume.
  //const double cubic_cell_volume=lattice_constant*lattice_constant*lattice_constant;
  const double primitive_cell_volume = lattice_vectors[0].cross(lattice_vectors[1]).dot(lattice_vectors[2]);

  double cFideal_of_primitive_cell=0;
  {
    //Find inhomogeneous Fideal of one crystal primitive cell
    // scale our dx by w.
    const double dx = dx_input*gwidth;
    printf("using dx = %g for ideal gas free energy integral\n", dx);

    const int Nl = (lattice_constant/2)/dx+1; // number of infinitesimal lengths along one of the lattice vectors
    printf("New length Nl*dx is off from a/2 by a factor of (Nl*dx)/(a/2)=%g\n", (Nl*dx)/(lattice_constant/2));  //Looking for source of errors
    //Nl^3 is total number of infinitesimal parallelepipeds (of volume dV) in one primitive cell
    //The volume of one infinitesimal parallelepiped dV=2dx^3 with the current definition of dx="dx_proper"/2
    //where dV=(dx_proper^3)/4 just as V=(a^3)/4 is the volume of one parallelepiped with lattice_constant=a.
    //NOTE: dx_proper chops up lattice constant a, whereas the current definition of dx chops up a/2 to create
    //chunks that correspond in number to infinitesimal lengths, or chunks, along the lattice_vectors.
    //dx_proper=lattice_constant/number of chunks along lattice vector
    //dx=(lattice_constant/2)/number of chunks along lattice vector
    // const double dV = uipow(lattice_constant/Nl,3)/4.0;

    const vector3d da1 = lattice_vectors[0]/Nl; //infinitesimal lattice vectors 1/Nl of lattice vector
    const vector3d da2 = lattice_vectors[1]/Nl;
    const vector3d da3 = lattice_vectors[2]/Nl;
    printf("da1x = %g vs dx = %g\n", da2.x, dx);
    const double dV = da1.cross(da2).dot(da3); //volume of infinitesimal parallelpiped

    const double analytic_ideal_free_energy =
      (1-fv)*temp*(log((1-fv)*2.646476976618268e-6/sqrt(temp*temp*temp))
                   - 3*log(sqrt(2*M_PI)*gwidth)
                   - 5.0/2);
    if (gwidth < 0.01*lattice_constant) {
      printf("gwidth is very small, so I'm trusting our analytic ideal free energy.\n");
      cFideal_of_primitive_cell = analytic_ideal_free_energy;
      printf("analytic crystal ideal gas free energy per volume = %.12g\n",
             cFideal_of_primitive_cell/primitive_cell_volume);
    } else {

      for (int i=0; i<Nl; i++) {  //integrate over one primitive cell
        for (int j=0; j<Nl; j++) {
          for (int k=0; k<Nl; k++) {
            vector3d r=i*da1 + j*da2 + k*da3;

            const int many_cells=2 + 6*gwidth*sqrt(2)/lattice_constant; // how many cells to sum over   // original MS
            //const int many_cells=4 + 12*gwidth*sqrt(2)/lattice_constant; // how many cells to sum over   // testing higher temps 2021
            //Gaussians father away won't contriubute much
            double n = 0;
            const double kT = temp;
            for (int t=-many_cells; t <=many_cells; t++) {
              for(int u=-many_cells; u<=many_cells; u++)  {
                for (int v=-many_cells; v<= many_cells; v++) {
                  const vector3d R = t*lattice_vectors[0] + u*lattice_vectors[1] + v*lattice_vectors[2];
                  double deltar = (r-R).norm();
                  n += (1-fv)*exp(-deltar*deltar*(0.5/(gwidth*gwidth)))/uipow(sqrt(2*M_PI)*gwidth,3);
                }
              }
            }
            if (n > 1e-200) { // Only use n values that are large enough - avoid underflow and n=0 issues ln(0)=ERROR
              // printf("n = %g  dF = %g\n", n, dF);
              cFideal_of_primitive_cell += kT*n*(log(n*2.646476976618268e-6/(sqrt(kT)*kT)) - 1.0)*dV;
            }
          }
        }
      }  //End inhomogeneous Fideal calculation
      printf("crystal ideal gas free energy per volume = %.12g\n",
             cFideal_of_primitive_cell/primitive_cell_volume);
      printf("analytic ideal gas free energy per vol   = %g\n",
             analytic_ideal_free_energy/primitive_cell_volume);
    }
  }

  //Integrate over one primitive cell (a parallelepiped) to find free energy
  double cfree_energy_per_atom;
  double cfree_energy_per_vol;
  double hfree_energy_per_atom;
  double hfree_energy_per_vol;

  printf("\nCalculating Homogeneous Free Energy analytically ...\n");
  HomogeneousSFMTFluid hf = sfmt_homogeneous(reduced_density, temp);   //note: homogeneousFE/atom does not depend on fv or gw
  //Note: hf.energy() returns energy/volume

  hfree_energy_per_atom = (hf.energy()*primitive_cell_volume)/reduced_num_spheres;
  hfree_energy_per_vol = hf.energy();    // hf.energy() is free energy per vol
  hf.printme("     homogeneous:");
  printf("homogeneous free_energy per vol is %g\n", hf.energy());
  printf("homogeneous free_energy per atom is %g\n", hfree_energy_per_atom);

  // scale our dx by Xi or w, whichever is larger.
  const double dx = dx_input*(find_Xi(temp) + gwidth);
  printf("using dx = %g for excess free energy integral\n", dx);

  const int Nl = (lattice_constant/2)/dx+1; // number of infinitesimal lengths along one of the latt

  const vector3d da1 = lattice_vectors[0]/Nl; //infinitesimal lattice vectors 1/Nl of lattice vector
  const vector3d da2 = lattice_vectors[1]/Nl;
  const vector3d da3 = lattice_vectors[2]/Nl;
  printf("da1x = %g vs dx = %g\n", da2.x, dx);
  const double dV = da1.cross(da2).dot(da3); //volume of infinitesimal parallelpiped

  double cFexcess_of_primitive_cell=0;  //Crystal Excess Free Energy over one primitive cell
  double total_phi_1 = 0, total_phi_2 = 0, total_phi_3 = 0;

  const double max_distance_considered = radius_of_peak(gwidth, temp);
  const int many_cells = 2*max_distance_considered/lattice_constant+1;
  printf("XXXXXXXXXXXX many_cells is %d based on %g vs %g XXXXXXXXXXXXXX\n",
         many_cells, max_distance_considered, lattice_constant);
  printf("Nl is %d\n", Nl);

  const double alpha = find_alpha(temp);
  const double Xi = find_Xi(temp);

  double mean_n0 = 0, mean_n1 = 0, mean_n2 = 0, mean_n3 = 0;
  for (int i=0; i<Nl; i++) {
    for (int j=0; j<Nl; j++) {
      for (int k=0; k<Nl; k++) {
        vector3d r= i*da1 + j*da2 + k*da3;

        double n_0=0, n_1=0, n_2=0, n_3=0;  //weighted densities  (fundamental measures)
        vector3d nv_1, nv_2;
        tensor3d nm_2;
        weight n_weight;

        for (int t=-many_cells; t <=many_cells; t++) {
          for(int u=-many_cells; u<=many_cells; u++)  {
            for (int v=-many_cells; v<= many_cells; v++) {
              const vector3d R = t*lattice_vectors[0] + u*lattice_vectors[1] + v*lattice_vectors[2];
              if ((R-r).norm() < max_distance_considered) {
                if (MC_ERROR == 0) {
                  n_weight=find_weighted_den_aboutR_mc(r, R, dx, temp,
                                                       lattice_constant, gwidth, fv);
                } else {
                  if (my_experiment > 0) {
                    n_weight=find_weighted_den_aboutR_mc_accurately_experiment(r, R, gwidth, fv, alpha, Xi);
                  } else {
                    n_weight=find_weighted_den_aboutR_mc_accurately(r, R, gwidth, fv, alpha, Xi);
                  }
                  //double n3_error=report_my_error(r, R, gwidth, fv, alpha, Xi);  //FOR DEBUG - delete!
                  //printf(">>>>n3_error=%g\n",n3_error);   //FOR DEBUG - delete!
                  //long total_num_points=report_total_num_points(r, R, gwidth, fv, alpha, Xi);  //FOR DEBUG - delete!
                  //printf(">>>>total_num_points=%ld\n", total_num_points);   //FOR DEBUG - delete!
                }

                n_0 +=n_weight.n_0;
                n_1 +=n_weight.n_1;
                n_2 +=n_weight.n_2;
                n_3 +=n_weight.n_3;

                nv_1 +=n_weight.nv_1;
                nv_2 +=n_weight.nv_2;
                nm_2 += n_weight.nm_2;
              }
            }
          }
        }

        double phi_1 = -n_0*1.0*log(1.0-1.0*n_3);
        double phi_2 = (n_1*n_2 - nv_1.dot(nv_2))/(1-n_3);
        double phi_3;
        if (use_tensor_weight) {
          //double traceof_nm_2squared =  nm_2.x.x*nm_2.x.x +
                                        //nm_2.y.y*nm_2.y.y +
                                        //nm_2.z.z*nm_2.z.z +
                                      //2*nm_2.x.y*nm_2.y.x +
                                      //2*nm_2.x.z*nm_2.z.x +
                                      //2*nm_2.y.z*nm_2.z.y;

          double traceof_nm_2cubed =
            nm_2.x.x*nm_2.x.x*nm_2.x.x +   // xx xx xx
            nm_2.y.y*nm_2.y.y*nm_2.y.y +   // yy yy yy
            nm_2.z.z*nm_2.z.z*nm_2.z.z +   // zz zz zz
            3*nm_2.x.x*nm_2.x.y*nm_2.y.x + // xx xy yx three places to put xx
            3*nm_2.x.x*nm_2.x.z*nm_2.z.x + // xx xz zx three places to put xx
            3*nm_2.x.y*nm_2.y.y*nm_2.y.x + // xy yy yx three places to put yy
            3*nm_2.y.y*nm_2.y.z*nm_2.z.y + // yy yz zy three places to put yy
            3*nm_2.x.z*nm_2.z.z*nm_2.z.x + // xz zz zx three places to put zz
            3*nm_2.y.z*nm_2.z.z*nm_2.z.y + // yz zz zy three places to put zz
            3*nm_2.x.y*nm_2.y.z*nm_2.z.x + // xy yz zx three cyclic permutations
            3*nm_2.x.z*nm_2.z.y*nm_2.y.x;  // xz zy yx == the above, but these are the non-cyclic permutations

          //phi_3 = (uipow(n_2,3) - 3.0*n_2*nv_2.dot(nv_2) + (9/2.0)*(nv_2.dot(nm_2.dot(nv_2))-3.0*nm_2.determinant()))/(24*M_PI*uipow(1.0-1.0*n_3,2));        // Schmidt equation 5  (2000) - ok
          phi_3 =( n_2*(uipow(n_2,2) - 3.0*nv_2.dot(nv_2)) + (9/2.0)*(nv_2.dot(nm_2.dot(nv_2)) - traceof_nm_2cubed))
            /(24*M_PI*uipow(1-n_3,2));         // Santos (2012) - yes
          //phi_3 = (nv_2.dot(nm_2.dot(nv_2))-n_2*nv_2.dot(nv_2) - traceof_nm_2cubed + n_2*traceof_nm_2squared)/((16/3.0)*M_PI*uipow(1.0-1.0*n_3,2));      // Sweatman, Groh & Schmidt  (2001) - same as Santos with different w_nm
        } else {
          // The following was Rosenfelds early vector version of the functional
          //double phi_3 = (uipow(n_2,3) - 3*n_2*nv_2.normsquared())/(24*M_PI*uipow(1-n_3,2));
          // This is the fixed version, which comes from dimensional crossover
          phi_3 = uipow(n_2,3)*uipow(1.0 - nv_2.normsquared()/sqr(n_2),3)/(24*M_PI*sqr(1-n_3));
        }
        if (n_2 == 0 || sqr(n_2) <= nv_2.normsquared()) phi_3 = 0;
        if (n_0 < 1e-5*reduced_density && n_3 >= 1 && n_3 < 1 + 1e-14) {
          // This is a case where the prefactor on each of our free
          // energy contributions is within roundoff error of zero.
          // Sadly, n_3 may be greater than 1, in which case we would
          // get a NaN when the true contribution to the free energy
          // is probably actually zero.
          phi_1 = 0;
          phi_2 = 0;
          phi_3 = 0;
        }

        mean_n0 += n_0*dV;
        mean_n1 += n_1*dV;
        mean_n2 += n_2*dV;
        mean_n3 += n_3*dV;
        total_phi_1 += temp*phi_1*dV;
        total_phi_2 += temp*phi_2*dV;
        total_phi_3 += temp*phi_3*dV;
        cFexcess_of_primitive_cell += temp*(phi_1 + phi_2 + phi_3)*dV;  //NOTE: temp is Bolatzman constant times temperature divided by epsilon (which is 1)
        if (isnan(cFexcess_of_primitive_cell)) {                        // temp is a reduced temperature given by t* in the paper
          printf("free energy is a NaN!\n");
          printf("position is: %g %g %g\n", r.x, r.y, r.z);
          printf("n0 = %g\nn1 = %g\nn2=%g\nn3=%.17g = 1+%g\n", n_0, n_1, n_2, n_3, n_3-1);
          printf("phi1 = %g\nphi2 = %g\n\n#phi3=%g\n", phi_1, phi_2, phi_3);
          printf("             #NAN!\n"); //
          if (free_energy_output_file) {
            FILE *f = fopen(free_energy_output_file, "a");
            fprintf(f, "\n# git  version: %s\n", version_identifier());
            fprintf(f, "#dx\tmc-error\tphi1\tphi2\tphi3\tFtot\tFideal\t\t\tFEdiff\t\tmean_n3\t\t\trelerr\t\t\tmc_con\tmc_pf\tgw\t\tfv\tkT\tn\tseed\tmin\n");
            fprintf(f, "%g\t%g\t\t%g\t%g\t%g\t%g\t%g\t\t%g\t\t%g\t\t%g\t\t%ld\t%ld\t%g\t%g\t%g\t%g\t%g\t%g\n",
                    dx_input,
                    MC_ERROR,
                    total_phi_1/primitive_cell_volume,
                    total_phi_2/primitive_cell_volume,
                    total_phi_3/primitive_cell_volume,
                    cfree_energy_per_vol,
                    cFideal_of_primitive_cell/primitive_cell_volume,
                    total_phi_1+total_phi_2+total_phi_3,
                    mean_n3,
                    0.0,
                    mc_constant,
                    mc_prefactor,
                    gwidth,
                    fv,
                    temp,
                    reduced_density,
                    seed,
                    0.0);
            fclose(f);
          }
          data data_out;
          // Note that I use cFexcess_of_primitive_cell here which is
          // a NaN to represent the values that are actually
          // undefined.
          data_out.diff_free_energy_per_atom=cFexcess_of_primitive_cell;
          data_out.cfree_energy_per_atom=cFexcess_of_primitive_cell;
          data_out.hfree_energy_per_vol=hfree_energy_per_vol;
          data_out.cfree_energy_per_vol=cFexcess_of_primitive_cell;
          return data_out;
        }
        //printf("cFexcess_of_primitive_cell is now... %g\n", cFexcess_of_primitive_cell);   //debug
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
    //printf("cFexcess_of_primitive_cell so far=%g, phi_1=%g, phi_2=%g, phi_3=%g\n",cFexcess_of_primitive_cell, phi_1, phi_2, phi_3);
  }
  mean_n0 /= primitive_cell_volume;
  mean_n1 /= primitive_cell_volume;
  mean_n2 /= primitive_cell_volume;
  mean_n3 /= primitive_cell_volume;
  printf("mean n0 = %g, mean n1 = %g, mean n2 = %g, mean n3 = %g\n",
         mean_n0, mean_n1, mean_n2, mean_n3);
  printf("homo n0 = %g, homo n1 = %g, homo n2 = %g, homo n3 = %g\n",
         hf.get_n0(), hf.get_n1(), hf.get_n2(), hf.get_n3());

  double diff_n0=mean_n0- hf.get_n0();
  double diff_n1=mean_n1- hf.get_n1();
  double diff_n2=mean_n2- hf.get_n2();
  double diff_n3=mean_n3- hf.get_n3();
  printf("mean_n0-homo_n0 ratio=%g\n", diff_n0/hf.get_n0());    // HERE! Relative error
  printf("mean_n1-homo_n1 ratio=%g\n", diff_n1/hf.get_n1());
  printf("mean_n2-homo_n2 ratio=%g\n", diff_n2/hf.get_n2());
  printf("mean_n3-homo_n3 ratio=%g\n", diff_n3/hf.get_n3());

  // double array[4] = {diff_n0, diff_n1, diff_n2, diff_n3};
  // int i, j;
  // double tem;
  //  for (i=0; i < 4; ++i)
  //    if (array[i] < 0)
  //        array[i]=-array[i];  //take absolute value of difference
  //  //for (i=0; i < 4; ++i)
  //  //  printf("array[%i]: %g\n", i, array[i]);

  //  for (i=0; i < 4-1; ++i)   //sort to find greatest difference
  //    for (j=i+1; j<4; ++j)
  //       if (array[i] < array[j]) {
  //         tem = array[i];
  //         array[i] = array[j];
  //         array[j] = tem;
  //       }
  //  //for (i=0; i < 4; ++i)
  //  //  printf("array[%i]: %g\n", i, array[i]);

  // printf(">>>LARGEST diff in mean_n-homo_n=%g\n\n",array[0]);

  //printf("average of all mean_n-homo_n=%g\n", (mean_n0-hf.get_n0() + mean_n1-hf.get_n1()+mean_n2-hf.get_n2()+mean_n3-hf.get_n3())/4);// HERE! DELETE

  //There are 4 parallelepipeds in 1 cube; 1 atom/parallelepiped, 4 atoms/cube;
  //4*Vol_parallelepiped=Vol_cube=lattice_constant^3
  //cFexcess_of_primitive_cell and Fideal are over one parallelepiped (one primitive cell) with 1-fv atoms
  cfree_energy_per_atom=(cFideal_of_primitive_cell + cFexcess_of_primitive_cell)/reduced_num_spheres; //Fideal is the total inhomogeneous ideal free energy for 1 primitive cell
  cfree_energy_per_vol=(cFideal_of_primitive_cell + cFexcess_of_primitive_cell)/primitive_cell_volume; //
  printf("primitive cell volume = %g\n", primitive_cell_volume); //
  printf("cubic cell volume = %g;   cubic cell volume/4= %g\n", lattice_constant*lattice_constant*lattice_constant, (lattice_constant*lattice_constant*lattice_constant)/4); //
  printf("             phi_1 per volume = %g\n", total_phi_1/primitive_cell_volume); //
  printf("             phi_2 per volume = %g\n", total_phi_2/primitive_cell_volume); //
  printf("             phi_3 per volume = %g\n", total_phi_3/primitive_cell_volume); //

  printf("Crystal Ideal free energy per volume = %g\n", cFideal_of_primitive_cell/primitive_cell_volume);
  printf("Crystal Excess free energy per volume = %g\n", cFexcess_of_primitive_cell/primitive_cell_volume);
  printf("     Total crystal free energy per volume = %g\n", cfree_energy_per_vol);
  printf("cFideal_of_primitive_cell = %g\n", cFideal_of_primitive_cell);

  data data_out;
  data_out.diff_free_energy_per_atom=cfree_energy_per_atom-hfree_energy_per_atom;
  data_out.cfree_energy_per_atom=cfree_energy_per_atom;
  data_out.hfree_energy_per_vol=hfree_energy_per_vol;
  data_out.cfree_energy_per_vol=cfree_energy_per_vol;

  printf("*Homogeneous free energy calculated analytically\n");

  printf("data_out is: homFEperatom=%g, cryFEperatom=%g\n", hfree_energy_per_atom, data_out.cfree_energy_per_atom);
  printf("data_out is: homFEpervol=%g, cryFEpervol=%g\n", data_out.hfree_energy_per_vol, data_out.cfree_energy_per_vol);
  printf("data_out is: diffperatom=%g\n", data_out.diff_free_energy_per_atom);

  int create_alldat_file = 1;  //set to 0 for no alldat file, set to 1 to create alldat file
  // that saves data for every combination of gw and fv values ran
  double run_time = time() - start_time; // run time in seconds
  if (create_alldat_file > 0)  {
    // Create all output data filename
    char *alldat_filename = new char[1024];
    char *alldat_filedescriptor = new char[1024];
    sprintf(alldat_filedescriptor, "kT%5.3f_n%05.3f_fv%04.2f_gw%04.3f",     //original  MS
    //sprintf(alldat_filedescriptor, "kT%5.3f_n%05.3f_fv%04.2f_gw%06.5f",       //testing higher temps 2021
            temp, reduced_density, fv, gwidth);
    if (use_tensor_weight)  {
      sprintf(alldat_filename, "%s/%s-alldat_tensor.dat", data_dir, alldat_filedescriptor);
    } else {
      sprintf(alldat_filename, "%s/%s-alldat.dat", data_dir, alldat_filedescriptor);
    }
    printf("Create data file: %s\n", alldat_filename);

    //Create dataout file
    FILE *newmeltoutfile = fopen(alldat_filename, "w");
    if (newmeltoutfile) {
      fprintf(newmeltoutfile, "# git  version: %s\n", version_identifier());
      fprintf(newmeltoutfile, "#kT\tn\tfv\tgwidth\thFE/atom\tcFE/atom\tFEdiff/atom\tlat_const\tNsph\tdx\tmcconstant\tmcprefactor\tmcerror\tmcseed\ttime(h)\ttensor\n");
      fprintf(newmeltoutfile, "%g\t%g\t%g\t%g\t%g\t%g\t      %g\t\t%g\t\t%g\t%g\t%li\t%li\t%g\t%g\t%g\t%d\n",
              temp, reduced_density, fv, gwidth, hfree_energy_per_atom,
              cfree_energy_per_atom, data_out.diff_free_energy_per_atom,
              lattice_constant, reduced_num_spheres, dx_input, mc_constant, mc_prefactor, MC_ERROR, seed,
              run_time/60/60, use_tensor_weight);
      fclose(newmeltoutfile);
    } else {
      printf("Unable to open file %s!\n", alldat_filename);
    }
    delete[] alldat_filename;
  }
  printf("run time is %g hours\n", run_time/60/60);
  //printf("mc_constant=%ld, mc_prefactor=%ld\n", mc_constant, mc_prefactor); // HERE!

  if (free_energy_output_file) {
    FILE *f = fopen(free_energy_output_file, "a");
    fprintf(f, "\n# git  version: %s\n", version_identifier());
    fprintf(f, "#dx\tmc-error\tphi1\tphi2\tphi3\tFtot\tFideal\t\t\tFEdiff\t\tmean_n3\t\t\trelerr\t\t\tmc_con\tmc_pf\tgw\t\tfv\tkT\tn\tseed\tmin\ttensor\n");
    fprintf(f, "%g\t%g\t\t%g\t%g\t%g\t%g\t%g\t\t%g\t\t%g\t\t%g\t\t%ld\t%ld\t%g\t%g\t%g\t%g\t%g\t%g\t%d\n",
            dx_input,
            MC_ERROR,
            total_phi_1/primitive_cell_volume,
            total_phi_2/primitive_cell_volume,
            total_phi_3/primitive_cell_volume,
            cfree_energy_per_vol,
            cFideal_of_primitive_cell/primitive_cell_volume,
            data_out.diff_free_energy_per_atom,
            mean_n3,
            diff_n3/hf.get_n3(),
            mc_constant,
            mc_prefactor,
            gwidth,
            fv,
            temp,
            reduced_density,
            seed,
            run_time/60,
            use_tensor_weight);
    fclose(f);
  }
  if (my_experiment > 0) {
    printf("My experiment\n");
  }
  //printf("scaled num_points=%g\n", 20+200*gwidth); // HERE!
  //printf("scaled num_points=50000\n"); // HERE!
  return data_out;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%END NEW ENERGY FUNCTION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    data dhill_data=find_energy_new(temp, reduced_density, simplex_fe[k][0], simplex_fe[k][1], data_dir, dx, verbose);
    if (isnan(dhill_data.cfree_energy_per_vol)) {   //If crystal energy is NaN
      simplex_fe[k][2]=1e100;
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
  //double seed=1;
  double reduced_density=1.0, gw=-1, fv=-1, temp=1.0;
  //gw is standard deviation of a Gaussian function typically called "sigma" (not to be confused with the sigma of the WCA potential!)
  //in this program the term sigma means the sigma of the WCA potential
  //reduced_density is the homogeneous (flat) number density accounting for sphere vacancies MULTIPLIED BY the WCA sigma^3 to make it dimensionless
  //temp is the Boltzman constant MULTIPLIED BY the temperature in Kelvin


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
  sprintf(data_dir,"crystallization");


  //********************Setup POPT to get inputs from command line*******************

  // ----------------------------------------------------------------------------
  // Parse input options
  // ----------------------------------------------------------------------------

  bool show_version = false;
  poptOption optionsTable[] = {

    {
      "verbose", '\0', POPT_ARG_NONE, &verbose, 0,
      "Print lots of good stuff!", "BOOLEAN"
    },
    {
      "version", '\0', POPT_ARG_NONE, &show_version, 0,
      "Show the git-describe version", "BOOLEAN"
    },

    /*** FLUID PARAMETERS ***/
    {"kT", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &temp, 0, "temperature", "DOUBLE"},
    {"n", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &reduced_density, 0, "number density times WCA sigma^3", "DOUBLE"},
    {"fv", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &fv, 0, "fraction of vacancies", "DOUBLE"},
    {"gw", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &gw, 0, "width of Gaussian", "DOUBLE"},

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
    {
      "dx", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &dx, 0,
      "dimensionless grid spacing dx (scaled by w or Xi)", "DOUBLE"
    },

    /*** MONTE-CARLO SEED OPTIONS ***/
    {"mc", '\0', POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT, &NUM_POINTS, 0, "Number of Points for Monte-Carlo", "INT"},
    {"mc-error", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &MC_ERROR, 0, "Desired error in Monte-Carlo", "DOUBLE"},
    {"seed", '\0', POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT, &seed, 0, "Monte-Carlo seed", "DOUBLE"},

    {"mc-prefactor", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &mc_prefactor, 0, "Monte-Carlo mc-prefactor", "LONG"},   //temporary - delete later!
    {"mc-constant", '\0', POPT_ARG_LONG | POPT_ARGFLAG_SHOW_DEFAULT, &mc_constant, 0, "Monte-Carlo mc-constant", "LONG"},   //temporary - delete later!

    /*** PARAMETERS DETERMINING OUTPUT FILE DIRECTORY AND NAMES ***/
    // Where to save the free energy info
    {
      "filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &free_energy_output_file, 0,
      "File to append free energy info to", "FNAME"
    },
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
  if (show_version) {
    printf("new-melting version %s (should match `git describe --dirty`)\n", version_identifier());
    exit(0);
  }

  random::seed(seed);

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
  printf("\nUsing data directory: [deft/papers/fuzzy-fmt]/%s\n", data_dir);
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

    printf("\n\nreduced_density = %g, fv = %g, gw = %g  alpha=%g Xi=%g\n", reduced_density, fv, gw,
           find_alpha(temp), find_Xi(temp));

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
    double Xi=find_Xi(temp);
    double norm=uipow(gw*sqrt(2*M_PI),3);
    printf("Analytic Solution:\n");
    //n_0
    double n_0_of_r_1stterm = (uipow(gw*sqrt(2*M_PI),3)/(4*M_PI*r.z*r.z*norm*Xi*sqrt(M_PI)))*exp(-uipow((r.z-(alpha/2))/Xi,2));   //Uses first Term of w1 Taylor expansion
    double n_0_of_r = n_0_of_r_1stterm;
    printf("n_0   analytic n_0 = %g  for r.z=%g (compare with quadrature %g or mc %g)\n",
           n_0_of_r, r.z, w_R.n_0, w_MC.n_0);
    //printf("      1st term of n_0= %g\n", n_0_of_r_1stterm);
    //n_1
    double n_1_of_r_1stterm = (uipow(gw*sqrt(2*M_PI),3)/(4*M_PI*r.z*norm*Xi*sqrt(M_PI)))*exp(-uipow((r.z-(alpha/2))/Xi,2));   //Uses first Term of w1 Taylor expansion
    double n_1_of_r = n_1_of_r_1stterm;
    printf("n_1   analytic n_1 = %g  for r.z=%g (compare with quadrature %g or mc %g)\n",
           n_1_of_r, r.z, w_R.n_1, w_MC.n_1);
    //printf("      1st term of n_1= %g\n", n_1_of_r_1stterm);
    //n_2
    double n_2_of_r_1stterm = (uipow(gw*sqrt(2*M_PI),3)/(norm*Xi*sqrt(M_PI)))*exp(-uipow((r.z-(alpha/2))/Xi,2));   //Uses first Term of w2 Taylor expansion
    double n_2_of_r_3rdterm = 6*sqrt(2)*M_PI*gw*gw*gw*gw*gw/(norm*Xi*Xi*r.z)*exp(-uipow((r.z-(alpha/2))/Xi,2))*(2*uipow((r.z-(alpha/2))/Xi,2)-1);
    double n_2_of_r = n_2_of_r_1stterm + n_2_of_r_3rdterm;   //2nd term (and all even terms) = zero
    printf("n_2   analytic n_2 = %g  for r.z=%g (compare with quadrature %g or mc %g)\n",
           n_2_of_r, r.z, w_R.n_2, w_MC.n_2);
    //printf("      1st term of n_2= %g   2nd term of n_2= 0   3rd term of n_2= %g\n", n_2_of_r_1stterm, n_2_of_r_3rdterm);
    //n_3
    double n_3_of_r_1stterm = (uipow(gw*sqrt(2*M_PI),3)*(1-erf((r.z-(alpha/2))/Xi))/(2*norm));   //Uses first Term of w1 Taylor expansion
    double n_3_of_r = n_3_of_r_1stterm;
    printf("n_3   analytic n_3 = %g  for r.z=%g (compare with quadrature %g or mc %g)\n",
           n_3_of_r, r.z, w_R.n_3, w_MC.n_3);
    //printf("      1st term of n_3= %g\n", n_3_of_r_1stterm);

    //printf("alpha = %g,  Xi=%g, temp=%g\n", alpha, Xi, temp);

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

  //printf("mc-error=%lg\n", MC_ERROR)
  //printf("mc-prefactor=%ld\n", mc_prefactor);
  //printf("mc-constant=%ld\n", mc_constant);

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
    if (use_tensor_weight)  {
      sprintf(bestdat_filename, "%s/kT%05.3f_n%05.3f_best_tensor.dat", data_dir, temp, reduced_density);
    } else {
      sprintf(bestdat_filename, "%s/kT%05.3f_n%05.3f_best.dat", data_dir, temp, reduced_density);
    }

    //Create bestdataout file
    printf("Create best data file: %s\n", bestdat_filename);
    FILE *newmeltbest = fopen(bestdat_filename, "w");
    if (newmeltbest) {
      fprintf(newmeltbest, "# git version: %s\n", version_identifier());

      fprintf(newmeltbest, "#kT\tn\tfv\tgwidth\thFE/atom\tbest_cFE/atom\tbest_FEdiff/atom\tbest_lat_const\tNsph\tdx\tmcerror\tmcseed\thFE/volume\tbest_cFE/volume\tmcconstant\tmcprefactor\ttensor\n");
      fprintf(newmeltbest, "%g\t%g\t%g\t%g\t%g\t%g\t\t%g\t\t%g\t\t%g\t%g\t%g\t%g\t%g\t%g\t%li\t%li\t%d\n", temp, reduced_density, best_fv, best_gwidth,
              best_cfree_energy-best_energy_diff, best_cfree_energy, best_energy_diff,
              best_lattice_constant, 1-fv, dx, MC_ERROR, seed, hfree_energy_pervol, cfree_energy_pervol, mc_constant, mc_prefactor, use_tensor_weight);    //Nsph=1-fv for parallepiped
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
    if (use_tensor_weight) {
      sprintf(bestdat_filename, "%s/kT%05.3f_n%05.3f_best_tensor.dat", data_dir, temp, reduced_density);
    } else {
      sprintf(bestdat_filename, "%s/kT%05.3f_n%05.3f_best.dat", data_dir, temp, reduced_density);
    }

    //Create bestdataout file
    printf("Create best data file: %s\n", bestdat_filename);
    FILE *newmeltbest = fopen(bestdat_filename, "w");
    if (newmeltbest) {
      fprintf(newmeltbest, "# git version: %s\n", version_identifier());
      fprintf(newmeltbest, "#kT\tn\tfv\tgwidth\thFE/atom\tbest_cFE/atom\tbest_FEdiff/atom\tbest_lat_const\tNsph\tdx\tmcerror\tmcseed\thFE/volume\tbest_cFE/volume\tmcconstant\tmcprefactor\ttensor\n");
      fprintf(newmeltbest, "%g\t%g\t%g\t%g\t%g\t%g\t\t%g\t\t%g\t\t%g\t%g\t%g\t%g\t%g\t%g\t%li\t%li\t%d\n", temp, reduced_density, best_fv, best_gwidth,
              best_cfree_energy-best_energy_diff, best_cfree_energy, best_energy_diff,
              best_lattice_constant, 1-fv, dx, MC_ERROR, seed, hfree_energy_pervol, cfree_energy_pervol, mc_constant, mc_prefactor, use_tensor_weight);    //Nsph=1-fv for parallepiped
      fclose(newmeltbest);
    } else {
      printf("Unable to open file %s!\n", bestdat_filename);
    }
    delete[] bestdat_filename;

  } else {
    find_energy_new(temp, reduced_density, fv, gw, data_dir, dx, bool(verbose));
  }

  return 0;
}
