// -*- mode: C++; -*-

#include "new/NewFunctional.h"
#include "utilities.h"
#include "handymath.h"

#pragma GCC diagnostic ignored "-Wunused-variable"
class SW_liquid : public NewFunctional {
public:
SW_liquid(int myNx, int myNy, int myNz);
SW_liquid(double ax, double ay, double az, double dx);
double true_energy() const;
Vector grad() const;
bool have_preconditioner() const;
EnergyGradAndPrecond energy_grad_and_precond() const;
void printme(const char *prefix) const;
double get_ESW() const;
double get_Fideal() const;
Vector get_dAdR() const;
double get_dV() const;
ComplexVector get_deltak() const;
ComplexVector get_deltak2() const;
ComplexVector get_deltaprime() const;
ComplexVector get_deltaprimex() const;
ComplexVector get_deltaprimey() const;
ComplexVector get_deltaprimez() const;
ComplexVector get_deltax() const;
ComplexVector get_deltay() const;
ComplexVector get_deltaz() const;
Vector get_dn1v_dot_n2v_by_dn1vx() const;
Vector get_dn1v_dot_n2v_by_dn1vy() const;
Vector get_dn1v_dot_n2v_by_dn1vz() const;
Vector get_dn1v_dot_n2v_by_dn2vx() const;
Vector get_dn1v_dot_n2v_by_dn2vy() const;
Vector get_dn1v_dot_n2v_by_dn2vz() const;
Vector get_dn2vsqr_by_dn2vx() const;
Vector get_dn2vsqr_by_dn2vy() const;
Vector get_dn2vsqr_by_dn2vz() const;
Vector get_dphi1_by_dn0() const;
Vector get_dphi1_by_dn3() const;
Vector get_dphi2_by_dn1() const;
Vector get_dphi2_by_dn1vx() const;
Vector get_dphi2_by_dn1vy() const;
Vector get_dphi2_by_dn1vz() const;
Vector get_dphi2_by_dn2() const;
Vector get_dphi2_by_dn2vx() const;
Vector get_dphi2_by_dn2vy() const;
Vector get_dphi2_by_dn2vz() const;
Vector get_dphi2_by_dn3() const;
Vector get_dphi3_by_dn2() const;
Vector get_dphi3_by_dn2vx() const;
Vector get_dphi3_by_dn2vy() const;
Vector get_dphi3_by_dn2vz() const;
Vector get_dphi3_by_dn3() const;
Vector get_dphitot_by_dn0() const;
Vector get_dphitot_by_dn1() const;
Vector get_dphitot_by_dn1vx() const;
Vector get_dphitot_by_dn1vy() const;
Vector get_dphitot_by_dn1vz() const;
Vector get_dphitot_by_dn2() const;
Vector get_dphitot_by_dn2vx() const;
Vector get_dphitot_by_dn2vy() const;
Vector get_dphitot_by_dn2vz() const;
Vector get_dphitot_by_dn3() const;
double get_dr() const;
double get_external() const;
Vector get_gSigmaA() const;
double get_kTphi1() const;
double get_kTphi2() const;
double get_kTphi3() const;
Vector get_n0() const;
Vector get_n1() const;
Vector get_n1v_dot_n2v() const;
Vector get_n1vx() const;
Vector get_n1vy() const;
Vector get_n1vz() const;
Vector get_n2() const;
Vector get_n2vsqr() const;
Vector get_n2vx() const;
Vector get_n2vy() const;
Vector get_n2vz() const;
Vector get_n3() const;
ComplexVector get_ngphi0() const;
ComplexVector get_ngphi1() const;
ComplexVector get_ngphi2() const;
ComplexVector get_ngphi3() const;
ComplexVector get_ngphi4() const;
Vector get_phi1() const;
Vector get_phi2() const;
Vector get_phi3() const;
ComplexVector get_step() const;
double get_sw() const;
double get_volume() const;
double get_whitebear() const;
ComplexVector get_xi0phik() const;
ComplexVector get_xi1phik() const;
ComplexVector get_xi2phik() const;
ComplexVector get_xi3phik() const;
ComplexVector get_xi4phik() const;
double d_by_dNx() const;
double d_by_dNy() const;
double d_by_dNz() const;
double d_by_dR() const;
double d_by_da1() const;
double d_by_da2() const;
double d_by_da3() const;
double d_by_depsilon() const;
double d_by_dkT() const;
double d_by_dlambda() const;
double d_by_dmu() const;
double d_by_dsigma() const;
Vector get_rx() const;
Vector get_ry() const;
Vector get_rz() const;
Vector get_r() const;
ComplexVector get_kx() const;
ComplexVector get_ky() const;
ComplexVector get_kz() const;
ComplexVector get_k() const;
double &Nx() const {
	int sofar = 0;
	return data[sofar];
}

double &Ny() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	return data[sofar];
}

double &Nz() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	return data[sofar];
}

double &R() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	return data[sofar];
}

double &a1() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	return data[sofar];
}

double &a2() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	return data[sofar];
}

double &a3() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	return data[sofar];
}

Vector Vext() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	return data.slice(sofar,Nx*Ny*Nz);
}

Vector n() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	sofar += Nx*Ny*Nz; // Vext
	return data.slice(sofar,Nx*Ny*Nz);
}

double &epsilon() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	sofar += Nx*Ny*Nz; // Vext
	sofar += Nx*Ny*Nz; // n
	return data[sofar];
}

double &kT() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	sofar += Nx*Ny*Nz; // Vext
	sofar += Nx*Ny*Nz; // n
	sofar += 1; // epsilon
	return data[sofar];
}

double &lambda() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	sofar += Nx*Ny*Nz; // Vext
	sofar += Nx*Ny*Nz; // n
	sofar += 1; // epsilon
	sofar += 1; // kT
	return data[sofar];
}

double &mu() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	sofar += Nx*Ny*Nz; // Vext
	sofar += Nx*Ny*Nz; // n
	sofar += 1; // epsilon
	sofar += 1; // kT
	sofar += 1; // lambda
	return data[sofar];
}

double &sigma() const {
	int sofar = 0;
	const double Nx = data[sofar++];
	const double Ny = data[sofar++];
	const double Nz = data[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	sofar += Nx*Ny*Nz; // Vext
	sofar += Nx*Ny*Nz; // n
	sofar += 1; // epsilon
	sofar += 1; // kT
	sofar += 1; // lambda
	sofar += 1; // mu
	return data[sofar];
}

private:
	mutable double ESW;
	mutable double Fideal;
	mutable double dV;
	mutable double dr;
	mutable double external;
	mutable double kTphi1;
	mutable double kTphi2;
	mutable double kTphi3;
	mutable double sw;
	mutable double volume;
	mutable double whitebear;
}; // End of SW_liquid class
#pragma GCC diagnostic pop
