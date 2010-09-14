#include <Eigen/LU>
#include <Eigen/Geometry>
#include "lattice.h"

Lattice::Lattice(Cartesian a1, Cartesian a2, Cartesian a3) {
  R.col(0) = a1; R.col(1) = a2; R.col(2) = a3;
  Rinv = R.inverse();
  vol = R.determinant();

  // Find a "minimal" basis set, that cannot be overly-skewed.
  minimalR = R;
  minimalRinv = Rinv;
  {
    bool change;
    const double MIN_SYMM_TOL = 1e-10;
    // First work on Wigner-Seitz cell...
    do {
      change = false;
      for (int k1 = 0; k1 < 3; k1 ++) {
        int k2 = (k1 + 1)%3;
        int k3 = (k1 + 2)%3;
        for (int i = -1; i <= 1; i++)
          for (int j = -1; j <= 1; j++) {
            Matrix3d d;
            d.setZero();
            d(k1,k1) = 1.0;
            d(k2,k1) = i;
            d(k3,k1) = j;
            d(k2,k2) = 1.0;
            d(k3,k3) = 1.0;
            Matrix3d newR = minimalR * d;
            if ((newR.transpose()*newR).trace() + MIN_SYMM_TOL
                < (minimalR.transpose()*minimalR).trace()) {
              change = true;
              minimalR = newR;
              //printf("newR.determinant() is %g\n", newR.determinant());
            }
          }
      }
    } while (change);
    // Now work on Brillouin zone...
    minimalRinv = minimalR.inverse();
  }
  //printf("R.determinant() is %g\n", R.determinant());
  //printf("minimalR.determinant() is %g\n", minimalR.determinant());
  assert(fabs(vol - minimalR.determinant()) < 1e-10*vol);
  assert(volume() > 0);
}

void Lattice::reorientBasis(Cartesian zdir) {
  zdir = zdir / zdir.norm();
  Cartesian a1new(1e10, 0, 0.3e9),
    a2new(1.3e9, 1e10, 0.7e9),
    a3new(1.23e9, 1e9, 0.73e9);
  bool got1 = false, got2 = false;
  for (int i=-3; i<=3; i++) {
    for (int j=-3; j<=3; j++) {
      for (int k=-3; k<=3; k++) {
        Cartesian here(i*R.col(0) + j*R.col(1) + k*R.col(2));
        double sin1 = here.cross(a1new).norm()/(here.norm()*a1new.norm());
        double sin2 = here.cross(a2new).norm()/(here.norm()*a2new.norm());
        double cosz = here.dot(zdir)/here.norm();
        if (fabs(cosz) <= 1e-7) {
          if (!got1) {
            a1new = here;
            got1 = true;
          } else if (!got2 && sin1 >= 1e-7) {
            a2new = here;
            got2 = true;
          } else if (here.norm() < a1new.norm() && sin2 >= 1e-7) {
            a1new = here;
          } else if (here.norm() < a2new.norm() && sin1 >= 1e-7) {
            a2new = here;
          }
        } else {
          if (here.norm() > 0 &&
              fabs(here.dot(zdir)) < fabs(a3new.dot(zdir)))
            a3new = here;
        }
      }
    }
  }
  // set up the new lattice...
  double vol0 = volume();
  R.col(0) = a1new; R.col(1) = a2new; R.col(2) = a3new;
  if (R.determinant() < 0) R.col(2) = -a3new;
  Rinv = R.inverse();
  vol = R.determinant();
  assert(volume() > 0);
  assert(fabs(volume()-vol0)/vol0 < 1e-7);
}

Cartesian Lattice::wignerSeitz(Cartesian r) const {
  // move_to_wigner_seitz_cell now can be called with a lattice that is not
  // minimal, and with a r vector that is not in the lattice parallelogon.
  // Note that r is input in cartesian coordinates, not lattice coordinates.
  
  // convert r to coefficients of minimized_lat vectors
  Relative relr(minimalRinv*r);
  
  // shift r to having components in [0, 1), so it is in the 
  // minimized parallelogon (works even for negative input components!)
  for (int i = 0; i < 3; ++i) relr(i) -= floor(relr(i)); 
  
  // convert back to cartesian coords
  r = minimalR * relr;

  // For anyone wondering about the insane unrolling here, it's an
  // optimization introduced because this code was (at one time)
  // significantly slowing down the JDFT fluid calculations. -- DJR
  // Note: If more optim is needed, cache what is subtracted from r. -- DHF
  Cartesian newr = r;
  double r2 = r.squaredNorm(), newr2;
  {
    const Cartesian rtry(r - minimalR.col(0)); // Cartesian(1,0,0);
    if ((newr2 = rtry.squaredNorm()) < r2) {
      newr = rtry;
      r2 = newr2;
    }
  }
  {
    const Cartesian rtry(r - minimalR.col(1)); // Cartesian(0,1,0);
    if ((newr2 = rtry.squaredNorm()) < r2) {
      newr = rtry;
      r2 = newr2;
    }
  }
  {
    const Cartesian rtry(r - minimalR.col(2)); // Cartesian(0,0,1);
    if ((newr2 = rtry.squaredNorm()) < r2) {
      newr = rtry;
      r2 = newr2;
    }
  }
  {
    const Cartesian rtry(r - minimalR*Cartesian(1,1,0));
    if ((newr2 = rtry.squaredNorm()) < r2) {
      newr = rtry;
      r2 = newr2;
    }
  }
  {
    const Cartesian rtry(r - minimalR*Cartesian(0,1,1));
    if ((newr2 = rtry.squaredNorm()) < r2) {
      newr = rtry;
      r2 = newr2;
    }
  }
  {
    const Cartesian rtry(r - minimalR*Cartesian(1,0,1));
    if ((newr2 = rtry.squaredNorm()) < r2) {
      newr = rtry;
      r2 = newr2;
    }
  }
  {
    const Cartesian rtry(r - minimalR*Cartesian(1,1,1));
    if (rtry.squaredNorm() < r2) newr = rtry;
  }
  return newr;
} // end wignerSeitz()

Reciprocal Lattice::brillouinZone(Reciprocal k) const {
  // move_to_wigner_seitz_cell now can be called with a lattice that is not
  // minimal, and with a r vector that is not in the lattice parallelogon.
  // Note that r is input in cartesian coordinates, not lattice coordinates.
  
  // convert k to coefficients of minimized_lat vectors
  RelativeReciprocal relk(minimalR.transpose()*k);
  
  // shift r to having components in [0, 1), so it is in the 
  // minimized parallelogon (works even for negative input components!)
  for (int i = 0; i < 3; ++i) relk(i) -= floor(relk(i)); 
  
  // convert back to cartesian coords
  k = minimalRinv.transpose() * relk;

  // For anyone wondering about the insane unrolling here, it's an
  // optimization introduced because this code was (at one time)
  // significantly slowing down the JDFT fluid calculations. -- DJR
  // Note: If more optim is needed, cache what is subtracted from r. -- DHF
  Reciprocal newk = k;
  double k2 = k.squaredNorm(), newk2;
  {
    const Reciprocal ktry(k - minimalR.col(0)); // Reciprocal(1,0,0);
    if ((newk2 = ktry.squaredNorm()) < k2) {
      newk = ktry;
      k2 = newk2;
    }
  }
  {
    const Reciprocal ktry(k - minimalR.col(1)); // Reciprocal(0,1,0);
    if ((newk2 = ktry.squaredNorm()) < k2) {
      newk = ktry;
      k2 = newk2;
    }
  }
  {
    const Reciprocal ktry(k - minimalR.col(2)); // Reciprocal(0,0,1);
    if ((newk2 = ktry.squaredNorm()) < k2) {
      newk = ktry;
      k2 = newk2;
    }
  }
  {
    const Reciprocal ktry(k - minimalR*Reciprocal(1,1,0));
    if ((newk2 = ktry.squaredNorm()) < k2) {
      newk = ktry;
      k2 = newk2;
    }
  }
  {
    const Reciprocal ktry(k - minimalR*Reciprocal(0,1,1));
    if ((newk2 = ktry.squaredNorm()) < k2) {
      newk = ktry;
      k2 = newk2;
    }
  }
  {
    const Reciprocal ktry(k - minimalR*Reciprocal(1,0,1));
    if ((newk2 = ktry.squaredNorm()) < k2) {
      newk = ktry;
      k2 = newk2;
    }
  }
  {
    const Reciprocal ktry(k - minimalR*Reciprocal(1,1,1));
    if (ktry.squaredNorm() < k2) newk = ktry;
  }
  return newk;
} // end brillouinZone()
