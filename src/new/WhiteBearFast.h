// -*- mode: C++; -*-

#include "new/NewFunctional.h"
#include "utilities.h"
#include "handymath.h"


class WhiteBear : public NewFunctional {
public:
WhiteBear(double varR_arg) : varR(varR_arg)  {
}

	double energy(const Vector &xxx) const;
	double energy_per_volume(const Vector &xxx) const;
	double denergy_per_volume_dx(const Vector &xxx) const;
	Vector grad(const Vector &xxx) const;
	void printme(const char * prefix) const;
	Vector createInput(double Nx, double Ny, double Nz, double R, double a1, double a2, double a3, Vector x, double kT) const;
double & Nx(Vector xxx) const {
	int sofar = 0;
	return xxx[sofar];
}

double & Ny(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	return xxx[sofar];
}

double & Nz(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	return xxx[sofar];
}

double & R(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	return xxx[sofar];
}

double & a1(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	return xxx[sofar];
}

double & a2(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	return xxx[sofar];
}

double & a3(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	return xxx[sofar];
}

Vector x(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	return xxx.slice(sofar,Nx*Ny*Nz);
}

double & kT(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	sofar += Nx*Ny*Nz; // x
	return xxx[sofar];
}

private:
	double varR;
	mutable double dV;
	mutable double dr;
	mutable double kTphi1;
	mutable double kTphi2;
	mutable double kTphi3;
	mutable double volume;
	mutable double whitebear;
}; // End of WhiteBear class
