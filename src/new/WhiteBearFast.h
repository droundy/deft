// -*- mode: C++; -*-

#include "new/Functional.h"
#include "utilities.h"
#include "handymath.h"


class WhiteBear : public Functional {
public:
WhiteBear(double R_arg) : R(R_arg)  {
}

	double energy(const Vector &xxx) const;
	double energy_per_volume(const Vector &xxx) const;
	double denergy_per_volume_dx(const Vector &xxx) const;
	Vector grad(const Vector &xxx) const;
	void printme(const char * prefix) const;
	Vector createInput(double Nx, double Ny, double Nz, double R, double a1, double a2, double a3, Vector x, double kT) const;
void set_Nx(Vector xxx, double Nx) const {
	int sofar = 0;
	xxx[sofar] = Nx;
}

void set_Ny(Vector xxx, double Ny) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	xxx[sofar] = Ny;
}

void set_Nz(Vector xxx, double Nz) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	xxx[sofar] = Nz;
}

void set_R(Vector xxx, double R) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	xxx[sofar] = R;
}

void set_a1(Vector xxx, double a1) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	xxx[sofar] = a1;
}

void set_a2(Vector xxx, double a2) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	xxx[sofar] = a2;
}

void set_a3(Vector xxx, double a3) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	xxx[sofar] = a3;
}

void set_x(Vector xxx, Vector x) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	xxx.slice(sofar,Nx*Ny*Nz) = x;
}

void set_kT(Vector xxx, double kT) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	sofar += 1; // a3
	sofar += Nx*Ny*Nz; // x
	xxx[sofar] = kT;
}

double get_Nx(Vector xxx) const {
	int sofar = 0;
	return xxx[sofar];
}

double get_Ny(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	return xxx[sofar];
}

double get_Nz(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	return xxx[sofar];
}

double get_R(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	return xxx[sofar];
}

double get_a1(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	return xxx[sofar];
}

double get_a2(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	return xxx[sofar];
}

double get_a3(Vector xxx) const {
	int sofar = 0;
	const double Nx = xxx[sofar++];
	const double Ny = xxx[sofar++];
	const double Nz = xxx[sofar++];
	sofar += 1; // R
	sofar += 1; // a1
	sofar += 1; // a2
	return xxx[sofar];
}

Vector get_x(Vector xxx) const {
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

double get_kT(Vector xxx) const {
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
	double R;
	mutable double dV;
	mutable double dr;
	mutable double volume;
	mutable double whitebear;
}; // End of WhiteBear class
