// -*- mode: C++; -*-

#include "MinimalFunctionals.h"
#include "utilities.h"
#include "handymath.h"


class WaterSaft_type : public FunctionalInterface {
public:
WaterSaft_type(double R_arg, double epsilon_association_arg, double kappa_association_arg, double epsilon_dispersion_arg, double lambda_dispersion_arg, double length_scaling_arg, double mu_arg) : R(R_arg), epsilon_association(epsilon_association_arg), kappa_association(kappa_association_arg), epsilon_dispersion(epsilon_dispersion_arg), lambda_dispersion(lambda_dispersion_arg), length_scaling(length_scaling_arg), mu(mu_arg)  {
	have_integral = true;
	// TODO: code to evaluate Fourier transforms goes here
	oldkT = 0.0/0.0;  // initialize to NaN so we'll have to define transforms

}
~WaterSaft_type() {

}

bool I_have_analytic_grad() const {
	return false;}

double integral(const GridDescription &gd, double kT, const VectorXd &x) const {
	if (oldkT != kT) {
		oldkT = kT;

	}
	double output=0;
	VectorXcd ktemp0(gd.NxNyNzOver2);
	ktemp0 = fft(gd, x);

	VectorXcd ktemp1(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp1[0] = ktemp0[i]*(50.26548245743669*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp1[i] = ktemp0[i]*(25.132741228718345*R*sin(2.0*R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1);
	}

	VectorXd rtemp2(gd.NxNyNz);
	rtemp2 = ifft(gd, ktemp1);

	ktemp1.resize(0); // KSpace
	VectorXcd ktemp3(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp3[0] = ktemp0[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp3[i] = ktemp0[i]*exp(-4.0*R*R*lambda_dispersion*lambda_dispersion*length_scaling*length_scaling*(0.5*k_i[2]*k_i[2] + 0.5*k_i[1]*k_i[1] + 0.5*k_i[0]*k_i[0]));
	}

	VectorXd rtemp4(gd.NxNyNz);
	rtemp4 = ifft(gd, ktemp3);

	ktemp3.resize(0); // KSpace
	VectorXcd ktemp5(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp5[0] = ktemp0[i]*(50.26548245743669*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp5[i] = ktemp0[i]*(25.132741228718345*R*sin(2.0*R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1);
	}

	VectorXd rtemp6(gd.NxNyNz);
	rtemp6 = ifft(gd, ktemp5);

	ktemp5.resize(0); // KSpace
	VectorXcd ktemp7(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp7[0] = ktemp0[i]*(12.566370614359172*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp7[i] = ktemp0[i]*(12.566370614359172*R*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1);
	}

	VectorXd rtemp8(gd.NxNyNz);
	rtemp8 = ifft(gd, ktemp7);

	ktemp7.resize(0); // KSpace
	VectorXcd ktemp9(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp9[0] = ktemp0[i]*(4.188790204786391*R*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp9[i] = ktemp0[i]*(12.566370614359172*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*t2*cos(t2) + sin(t2))/(sqrt(t1)*t1));
	}

	VectorXd rtemp10(gd.NxNyNz);
	rtemp10 = ifft(gd, ktemp9);

	ktemp9.resize(0); // KSpace
	VectorXd rtemp11(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp11[i] = rtemp8[i]/(1 + -1.0*rtemp10[i]);
	}

	VectorXcd ktemp12(gd.NxNyNzOver2);
	ktemp12 = fft(gd, rtemp11);

	VectorXcd ktemp13(gd.NxNyNzOver2);
	ktemp13[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp13[i] = ktemp0[i]*(12.566370614359172*complex(0,1)*k_i[0]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1);
	}

	VectorXd rtemp14(gd.NxNyNz);
	rtemp14 = ifft(gd, ktemp13);

	ktemp13.resize(0); // KSpace
	VectorXd rtemp15(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp15[i] = -1.0*rtemp14[i]/(1 + -1.0*rtemp10[i]);
	}

	VectorXcd ktemp16(gd.NxNyNzOver2);
	ktemp16 = fft(gd, rtemp15);

	rtemp15.resize(0); // Realspace
	VectorXcd ktemp17(gd.NxNyNzOver2);
	ktemp17[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp17[i] = ktemp0[i]*(12.566370614359172*complex(0,1)*k_i[1]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1);
	}

	VectorXd rtemp18(gd.NxNyNz);
	rtemp18 = ifft(gd, ktemp17);

	ktemp17.resize(0); // KSpace
	VectorXd rtemp19(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp19[i] = -1.0*rtemp18[i]/(1 + -1.0*rtemp10[i]);
	}

	VectorXcd ktemp20(gd.NxNyNzOver2);
	ktemp20 = fft(gd, rtemp19);

	rtemp19.resize(0); // Realspace
	ktemp0[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp0[i] = ktemp0[i]*(12.566370614359172*complex(0,1)*k_i[2]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1);
	}

	VectorXd rtemp22(gd.NxNyNz);
	rtemp22 = ifft(gd, ktemp0);

	ktemp0.resize(0); // KSpace
	VectorXd rtemp23(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp23[i] = -1.0*rtemp22[i]/(1 + -1.0*rtemp10[i]);
	}

	VectorXcd ktemp24(gd.NxNyNzOver2);
	ktemp24 = fft(gd, rtemp23);

	rtemp23.resize(0); // Realspace
	VectorXd rtemp25(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp25[i] = -1.0*log(1 + -1.0*rtemp10[i]);
	}

	VectorXcd ktemp26(gd.NxNyNzOver2);
	ktemp26 = fft(gd, rtemp25);

	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp10[i];
		const double t2 = rtemp14[i]*rtemp14[i];
		const double t3 = rtemp18[i]*rtemp18[i];
		const double t4 = rtemp22[i]*rtemp22[i];
		const double t5 = rtemp8[i]*rtemp8[i];
		const double t6 = 1/t1;
		const double t7 = 1/rtemp10[i];
		rtemp25[i] = t7*t7*rtemp11[i]*(t5 + -3.0*t4 + -3.0*t3 + -3.0*t2)*(t6*((-1.768388256576615e-2*t7 + 1.768388256576615e-2*t6)*(log(t1)/(t6*t6) + rtemp10[i]) + 8.841941282883075e-3) + -8.841941282883075e-3*1 + 1.768388256576615e-2*rtemp25[i]) + t6*(t6*(7.957747154594767e-2*t5 + -7.957747154594767e-2*t4 + -7.957747154594767e-2*t3 + -7.957747154594767e-2*t2) + 7.957747154594767e-2*rtemp8[i]/R)/R;
	}

	VectorXcd ktemp28(gd.NxNyNzOver2);
	ktemp28 = fft(gd, rtemp25);

	rtemp25.resize(0); // Realspace
	VectorXd rtemp29(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp10[i];
		rtemp29[i] = ((t1*t1*log(t1) + rtemp10[i])*(2.6525823848649224e-2*rtemp8[i]*rtemp8[i] + -2.6525823848649224e-2*rtemp22[i]*rtemp22[i] + -2.6525823848649224e-2*rtemp18[i]*rtemp18[i] + -2.6525823848649224e-2*rtemp14[i]*rtemp14[i])/(t1*rtemp10[i]*rtemp10[i]) + 7.957747154594767e-2*rtemp8[i]/R)/t1;
	}

	VectorXcd ktemp30(gd.NxNyNzOver2);
	ktemp30 = fft(gd, rtemp29);

	rtemp29.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp30[0] = -1.0*ktemp30[i]*(-25.132741228718345*R) + ktemp28[i]*(12.566370614359172*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sin(R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t2 = exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t3 = 1/(sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t4 = t1*t3;
		ktemp30[i] = -1.0*ktemp30[i]*(12.566370614359172*t2*(-1.0*R*cos(R/t3) + -1.0*t4)) + ktemp28[i]*(12.566370614359172*R*t2*t4);
	}

	ktemp28.resize(0); // KSpace
	VectorXd rtemp32(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp10[i];
		rtemp32[i] = rtemp14[i]*(-5.305164769729845e-2*rtemp11[i]*(t1*t1*log(t1) + rtemp10[i])/(rtemp10[i]*rtemp10[i]) + -7.957747154594767e-2*1/R)/t1;
	}

	VectorXcd ktemp33(gd.NxNyNzOver2);
	ktemp33 = fft(gd, rtemp32);

	rtemp32.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp33[0] = ktemp30[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp33[i] = ktemp33[i]*(12.566370614359172*R*complex(0,1)*k_i[0]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1) + ktemp30[i];
	}

	ktemp30.resize(0); // KSpace
	VectorXd rtemp35(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp10[i];
		rtemp35[i] = rtemp18[i]*(-5.305164769729845e-2*rtemp11[i]*(t1*t1*log(t1) + rtemp10[i])/(rtemp10[i]*rtemp10[i]) + -7.957747154594767e-2*1/R)/t1;
	}

	VectorXcd ktemp36(gd.NxNyNzOver2);
	ktemp36 = fft(gd, rtemp35);

	rtemp35.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp36[0] = ktemp33[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp36[i] = ktemp36[i]*(12.566370614359172*R*complex(0,1)*k_i[1]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1) + ktemp33[i];
	}

	ktemp33.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp10[i];
		rtemp11[i] = rtemp22[i]*(-5.305164769729845e-2*rtemp11[i]*(t1*t1*log(t1) + rtemp10[i])/(rtemp10[i]*rtemp10[i]) + -7.957747154594767e-2*1/R)/t1;
	}

	VectorXcd ktemp39(gd.NxNyNzOver2);
	ktemp39 = fft(gd, rtemp11);

	rtemp11.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp39[0] = ktemp36[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp39[i] = ktemp39[i]*(12.566370614359172*R*complex(0,1)*k_i[2]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1) + ktemp36[i];
	}

	ktemp36.resize(0); // KSpace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp39[0] = (-7.957747154594767e-2*ktemp12[i]*(-25.132741228718345*R) + (ktemp26[i]*(-0.15915494309189535*(12.566370614359172*R*R)/R + -7.957747154594767e-2*(-25.132741228718345*R)) + -7.957747154594767e-2*ktemp12[i]*(12.566370614359172*R*R))/R)/R + ktemp39[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = 12.566370614359172*R*sin(R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]))*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]))/(sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t2 = exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t3 = -1.0*sin(R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]))/(sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0])) + R*cos(R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t4 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t5 = 12.566370614359172*t2*(-1.0*sin(R*sqrt(t4))/(sqrt(t4)) + -1.0*R*cos(R*sqrt(t4)));
		const complex t6 = sin(R*sqrt(t4));
		ktemp39[i] = (7.957747154594767e-2*ktemp24[i]*(12.566370614359172*R*t2*t6*complex(0,1)*k_i[2]/(sqrt(t4))) + 7.957747154594767e-2*ktemp20[i]*(12.566370614359172*R*t2*t6*complex(0,1)*k_i[1]/(sqrt(t4))) + 7.957747154594767e-2*ktemp16[i]*(12.566370614359172*R*t2*t6*complex(0,1)*k_i[0]/(sqrt(t4))) + -7.957747154594767e-2*t5*ktemp12[i] + (-7.957747154594767e-2*ktemp26[i]*(2.0*t1/R + t5) + 7.957747154594767e-2*ktemp24[i]*(12.566370614359172*t2*t3*complex(0,1)*k_i[2]/t4) + 7.957747154594767e-2*ktemp20[i]*(12.566370614359172*t2*t3*complex(0,1)*k_i[1]/t4) + 7.957747154594767e-2*ktemp16[i]*(12.566370614359172*t2*t3*complex(0,1)*k_i[0]/t4) + -7.957747154594767e-2*t1*ktemp12[i])/R)/R + ktemp39[i];
	}

	ktemp26.resize(0); // KSpace
	ktemp24.resize(0); // KSpace
	ktemp20.resize(0); // KSpace
	ktemp16.resize(0); // KSpace
	ktemp12.resize(0); // KSpace
	VectorXd rtemp42(gd.NxNyNz);
	rtemp42 = ifft(gd, ktemp39);

	ktemp39.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp42[i] = rtemp42[i]/rtemp6[i];
	}

	rtemp6.resize(0); // Realspace
	double 	s44 = 0;
	for (int i=0; i<gd.NxNyNz; i++) {
		s44 += gd.dvolume*kT*x[i]*(-6.283185307179586*R*R*(sqrt(0.15915494309189535*kappa_association*rtemp2[i]*(0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp4[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp4[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp4[i]*rtemp4[i]*(73.49635953848843*R*R*R*rtemp4[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp4[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*rtemp4[i]*(73.49635953848843*R*R*R*rtemp4[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT + rtemp42[i])*(-1.0*1 + exp(epsilon_association/kT))/(R*R) + 1) + -1.0)/(kappa_association*rtemp2[i]*(0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp4[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp4[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp4[i]*rtemp4[i]*(73.49635953848843*R*R*R*rtemp4[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp4[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*rtemp4[i]*(73.49635953848843*R*R*R*rtemp4[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT + rtemp42[i])*(-1.0*1 + exp(epsilon_association/kT))) + 0.5*1 + log(12.566370614359172*R*R*(sqrt(0.15915494309189535*kappa_association*rtemp2[i]*(0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp4[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp4[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp4[i]*rtemp4[i]*(73.49635953848843*R*R*R*rtemp4[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp4[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*rtemp4[i]*(73.49635953848843*R*R*R*rtemp4[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT + rtemp42[i])*(-1.0*1 + exp(epsilon_association/kT))/(R*R) + 1) + -1.0)/(kappa_association*rtemp2[i]*(0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp4[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp4[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp4[i]*rtemp4[i]*(73.49635953848843*R*R*R*rtemp4[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp4[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*rtemp4[i]*(73.49635953848843*R*R*R*rtemp4[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT + rtemp42[i])*(-1.0*1 + exp(epsilon_association/kT)))));
	}

	rtemp42.resize(0); // Realspace
	rtemp2.resize(0); // Realspace
	Fassoc = 4.0*s44;
	double 	s45 = 0;
	for (int i=0; i<gd.NxNyNz; i++) {
		s45 += R*R*R*epsilon_dispersion*gd.dvolume*rtemp4[i]*x[i]*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)/((-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1));
	}

	a1integrated = -4.1887902047863905*s45;
	double 	s46 = 0;
	for (int i=0; i<gd.NxNyNz; i++) {
		s46 += R*R*R*epsilon_dispersion*epsilon_dispersion*gd.dvolume*rtemp4[i]*x[i]*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp4[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp4[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)*((-4.1887902047863905*R*R*R*rtemp4[i] + 1)*(-4.1887902047863905*R*R*R*rtemp4[i] + 1))*((-4.1887902047863905*R*R*R*rtemp4[i] + 1)*(-4.1887902047863905*R*R*R*rtemp4[i] + 1))/(kT*(4.0*(4.1887902047863905*R*R*R*rtemp4[i])*(4.1887902047863905*R*R*R*rtemp4[i]) + 16.755160819145562*R*R*R*rtemp4[i] + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp4[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp4[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp4[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1));
	}

	rtemp4.resize(0); // Realspace
	a2integrated = 2.0943951023931953*s46;
	Fdisp = a2integrated + a1integrated;
	double 	s47 = 0;
	for (int i=0; i<gd.NxNyNz; i++) {
		s47 += gd.dvolume*kT*x[i]*(-1.0*1 + log(2.6464769766182683e-6*x[i]/(sqrt(kT)*kT)));
	}

	Fideal = s47;
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp18[i] = -7.957747154594767e-2*rtemp18[i]*rtemp18[i] + -7.957747154594767e-2*rtemp14[i]*rtemp14[i];
	}

	rtemp14.resize(0); // Realspace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp18[i] = -7.957747154594767e-2*rtemp22[i]*rtemp22[i] + rtemp18[i];
	}

	rtemp22.resize(0); // Realspace
	double 	s50 = 0;
	for (int i=0; i<gd.NxNyNz; i++) {
		s50 += gd.dvolume*kT*(8.841941282883075e-3*rtemp8[i]*(log(1 + -1.0*rtemp10[i])*(1 + -1.0*rtemp10[i])*(1 + -1.0*rtemp10[i]) + rtemp10[i])*(rtemp8[i]*rtemp8[i] + 37.69911184307752*rtemp18[i])/(rtemp10[i]*rtemp10[i]*(1 + -1.0*rtemp10[i])*(1 + -1.0*rtemp10[i])) + ((7.957747154594767e-2*rtemp8[i]*rtemp8[i] + rtemp18[i])/(1 + -1.0*rtemp10[i]) + -7.957747154594767e-2*rtemp8[i]*log(1 + -1.0*rtemp10[i])/R)/R);
	}

	rtemp18.resize(0); // Realspace
	rtemp10.resize(0); // Realspace
	rtemp8.resize(0); // Realspace
	whitebear = s50;
	double 	s51 = 0;
	for (int i=0; i<gd.NxNyNz; i++) {
		s51 += gd.dvolume*mu*x[i];
	}

	FSAFT = whitebear + s51 + Fideal + Fdisp + Fassoc;
	output = FSAFT;
	// 20 Fourier transform used.
	// 17 temporaries made
	return output;

}

VectorXd transform(const GridDescription &gd, double kT, const VectorXd &x) const {
	assert(0);
}

double transform(double kT, double x) const {
	if (oldkT != kT) {
		oldkT = kT;

	}
	double output = 0;
		double 	n = x;
	double 	boltz = -1.0*1 + exp(epsilon_association/kT);
	double 	deltak = 12.566370614359172*R*R;
	double 	n1 = 7.957747154594767e-2*deltak*n/R;
	double 	step = 4.188790204786391*R*R*R;
	double 	n3 = n*step;
	double 	dphi2_by_dn2 = n1/(1 + -1.0*n3);
	double 	n2 = deltak*n;
	double 	deltax = 0;
	double 	n2vx = deltax*n;
	double 	deltay = 0;
	double 	n2vy = deltay*n;
	double 	deltaz = 0;
	double 	n2vz = deltaz*n;
	double 	n2vsqr = n2vz*n2vz + n2vy*n2vy + n2vx*n2vx;
	double 	dphi3_by_dn2 = (n2vsqr*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3)) + n2*n2*((log(1 + -1.0*n3)*(2.6525823848649224e-2*1/n3 + -5.305164769729845e-2) + 2.6525823848649224e-2)/n3 + 2.6525823848649224e-2*log(1 + -1.0*n3)))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphitot_by_dn2 = dphi3_by_dn2 + dphi2_by_dn2;
	double 	dn1v_dot_n2v_by_dn2vx = 7.957747154594767e-2*deltax*n/R;
	double 	dphi2_by_dn2vx = -1.0*dn1v_dot_n2v_by_dn2vx/(1 + -1.0*n3);
	double 	dn2vsqr_by_dn2vx = 2.0*n2vx;
	double 	dphi3_by_dn2vx = dn2vsqr_by_dn2vx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphitot_by_dn2vx = dphi3_by_dn2vx + dphi2_by_dn2vx;
	double 	dn1v_dot_n2v_by_dn2vy = 7.957747154594767e-2*deltay*n/R;
	double 	dphi2_by_dn2vy = -1.0*dn1v_dot_n2v_by_dn2vy/(1 + -1.0*n3);
	double 	dn2vsqr_by_dn2vy = 2.0*n2vy;
	double 	dphi3_by_dn2vy = dn2vsqr_by_dn2vy*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphitot_by_dn2vy = dphi3_by_dn2vy + dphi2_by_dn2vy;
	double 	dn1v_dot_n2v_by_dn2vz = 7.957747154594767e-2*deltaz*n/R;
	double 	dphi2_by_dn2vz = -1.0*dn1v_dot_n2v_by_dn2vz/(1 + -1.0*n3);
	double 	dn2vsqr_by_dn2vz = 2.0*n2vz;
	double 	dphi3_by_dn2vz = dn2vsqr_by_dn2vz*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphitot_by_dn2vz = dphi3_by_dn2vz + dphi2_by_dn2vz;
	double 	n0 = 7.957747154594767e-2*deltak*n/(R*R);
	double 	dphi1_by_dn3 = n0/(1 + -1.0*n3);
	double 	n1vx = 7.957747154594767e-2*deltax*n/R;
	double 	n1vy = 7.957747154594767e-2*deltay*n/R;
	double 	n1vz = 7.957747154594767e-2*deltaz*n/R;
	double 	n1v_dot_n2v = n1vz*n2vz + n1vy*n2vy + n1vx*n2vx;
	double 	dphi2_by_dn3 = (n1*n2 + -1.0*n1v_dot_n2v)/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphi3_by_dn3 = n2*(n2vsqr*(-5.305164769729845e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((-5.305164769729845e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(0.1061032953945969*1/(1 + -1.0*n3) + 5.305164769729845e-2))/(1 + -1.0*n3) + (2.6525823848649224e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(1 + -1.0*n3) + 5.305164769729845e-2) + 2.6525823848649224e-2)/n3)/n3) + n2*n2*(1.768388256576615e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((1.768388256576615e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(-3.53677651315323e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2))/(1 + -1.0*n3) + (-8.841941282883075e-3*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(1 + -1.0*n3) + -1.768388256576615e-2) + -8.841941282883075e-3)/n3)/n3))/(1 + -1.0*n3);
	double 	dphitot_by_dn3 = dphi3_by_dn3 + dphi2_by_dn3 + dphi1_by_dn3;
	double 	dphitot_by_dn0 = -1.0*log(1 + -1.0*n3);
	double 	dphitot_by_dn1 = n2/(1 + -1.0*n3);
	double 	dn1v_dot_n2v_by_dn1vx = n2vx;
	double 	dphitot_by_dn1vx = -1.0*dn1v_dot_n2v_by_dn1vx/(1 + -1.0*n3);
	double 	dn1v_dot_n2v_by_dn1vy = n2vy;
	double 	dphitot_by_dn1vy = -1.0*dn1v_dot_n2v_by_dn1vy/(1 + -1.0*n3);
	double 	dn1v_dot_n2v_by_dn1vz = n2vz;
	double 	dphitot_by_dn1vz = -1.0*dn1v_dot_n2v_by_dn1vz/(1 + -1.0*n3);
	double 	deltaprime = -25.132741228718345*R;
	double 	deltaprimex = 0;
	double 	deltaprimey = 0;
	double 	deltaprimez = 0;
	double 	dAdR = kT*n*(deltaprimez*dphitot_by_dn2vz + deltaprimey*dphitot_by_dn2vy + deltaprimex*dphitot_by_dn2vx + -1.0*deltaprime*dphitot_by_dn2 + deltak*dphitot_by_dn3 + (dphitot_by_dn1vz*(7.957747154594767e-2*deltaz/R + 7.957747154594767e-2*deltaprimez) + dphitot_by_dn1vy*(7.957747154594767e-2*deltay/R + 7.957747154594767e-2*deltaprimey) + dphitot_by_dn1vx*(7.957747154594767e-2*deltax/R + 7.957747154594767e-2*deltaprimex) + dphitot_by_dn1*(-7.957747154594767e-2*deltak/R + -7.957747154594767e-2*deltaprime) + dphitot_by_dn0*(-0.15915494309189535*deltak/R + -7.957747154594767e-2*deltaprime)/R)/R);
	double 	deltak2 = 50.26548245743669*R*R;
	double 	gSigmaA = dAdR/(deltak2*kT*n*n);
	double 	eta_d = 4.1887902047863905*R*R*R*n;
	double 	c1 = lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + 2.25855;
	double 	c2 = lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + -0.66927;
	double 	c3 = lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576;
	double 	eta_eff = eta_d*(eta_d*(c3*eta_d + c2) + c1);
	double 	ghs = (1 + -0.5*eta_eff)/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff));
	double 	deta_eff_by_deta_d = eta_d*(3.0*c3*eta_d + 2.0*c2) + c1;
	double 	dghs_by_deta_d = deta_eff_by_deta_d*((3.0*1 + -1.5*eta_eff)/(1 + -1.0*eta_eff) + -0.5)/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff));
	double 	da1_by_deta_d = epsilon_dispersion*(dghs_by_deta_d*eta_d + ghs)*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
	double 	dc1_by_dlambda_dispersion = -1.50349*1 + 0.498868*lambda_dispersion;
	double 	dc2_by_dlambda_dispersion = 1.40049*1 + -1.655478*lambda_dispersion;
	double 	dc3_by_dlambda_dispersion = -15.0427*1 + 10.61654*lambda_dispersion;
	double 	deta_eff_by_dlambda_dispersion = eta_d*(eta_d*(dc3_by_dlambda_dispersion*eta_d + dc2_by_dlambda_dispersion) + dc1_by_dlambda_dispersion);
	double 	dghs_by_dlambda_dispersion = deta_eff_by_dlambda_dispersion*((3.0*1 + -1.5*eta_eff)/(1 + -1.0*eta_eff) + -0.5)/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff));
	double 	da1_by_dlambda_dispersion = epsilon_dispersion*eta_d*(-12.0*ghs*lambda_dispersion*lambda_dispersion + dghs_by_dlambda_dispersion*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0));
	double 	gSW = (-8.333333333333333e-2*da1_by_dlambda_dispersion*lambda_dispersion/eta_d + 0.25*da1_by_deta_d)/kT + gSigmaA;
	double 	deltasaft = boltz*gSW*kappa_association;
	double 	delta2k = 50.26548245743669*R*R;
	double 	nA = 1.9894367886486918e-2*delta2k*n/(R*R);
	double 	X = (0.25*sqrt(8.0*deltasaft*nA + 1) + -0.25)/(deltasaft*nA);
	double 	Fassoc = kT*n*(2.0*1 + 4.0*log(X) + -2.0*X);
	double 	a1 = epsilon_dispersion*eta_d*ghs*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
	double 	a1integrated = a1*n;
	double 	KHS = (eta_d*(eta_d*(eta_d*(-4.0*1 + eta_d) + 6.0) + -4.0) + 1)/(eta_d*(4.0*1 + 4.0*eta_d) + 1);
	double 	a2 = 0.5*KHS*da1_by_deta_d*epsilon_dispersion*eta_d;
	double 	a2integrated = a2*n/kT;
	double 	Fdisp = a2integrated + a1integrated;
	double 	gpermol = 1822.8885;
	double 	mH2O = 18.01528*gpermol;
	double 	Fideal = kT*n*(-1.0*1 + log(15.74960994572242*n/(sqrt(kT)*kT*sqrt(mH2O)*mH2O)));
	double 	phi1 = -1.0*n0*log(1 + -1.0*n3);
	double 	phi2 = (n1*n2 + -1.0*n1v_dot_n2v)/(1 + -1.0*n3);
	double 	phi3 = n2*(n2vsqr*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3)) + n2*n2*((log(1 + -1.0*n3)*(8.841941282883075e-3*1/n3 + -1.768388256576615e-2) + 8.841941282883075e-3)/n3 + 8.841941282883075e-3*log(1 + -1.0*n3)))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	whitebear = kT*(phi3 + phi2 + phi1);
	double 	FSAFT = mu*n + whitebear + Fideal + Fdisp + Fassoc;
	output = FSAFT;

	return output;

}

bool append_to_name( std::string) const {
	return false;}

double derive(double kT, double x) const {
	if (oldkT != kT) {
		oldkT = kT;

	}
	double output = 0;
	double 	boltz = -1.0*1 + exp(epsilon_association/kT);
	double 	n = x;
	double 	deltak = 12.566370614359172*R*R;
	double 	n1 = 7.957747154594767e-2*deltak*n/R;
	double 	step = 4.188790204786391*R*R*R;
	double 	n3 = n*step;
	double 	dphi2_by_dn2 = n1/(1 + -1.0*n3);
	double 	n2 = deltak*n;
	double 	deltax = 0;
	double 	n2vx = deltax*n;
	double 	deltay = 0;
	double 	n2vy = deltay*n;
	double 	deltaz = 0;
	double 	n2vz = deltaz*n;
	double 	n2vsqr = n2vz*n2vz + n2vy*n2vy + n2vx*n2vx;
	double 	dphi3_by_dn2 = (n2vsqr*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3)) + n2*n2*((log(1 + -1.0*n3)*(2.6525823848649224e-2*1/n3 + -5.305164769729845e-2) + 2.6525823848649224e-2)/n3 + 2.6525823848649224e-2*log(1 + -1.0*n3)))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphitot_by_dn2 = dphi3_by_dn2 + dphi2_by_dn2;
	double 	dn1v_dot_n2v_by_dn2vx = 7.957747154594767e-2*deltax*n/R;
	double 	dphi2_by_dn2vx = -1.0*dn1v_dot_n2v_by_dn2vx/(1 + -1.0*n3);
	double 	dn2vsqr_by_dn2vx = 2.0*n2vx;
	double 	dphi3_by_dn2vx = dn2vsqr_by_dn2vx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphitot_by_dn2vx = dphi3_by_dn2vx + dphi2_by_dn2vx;
	double 	dn1v_dot_n2v_by_dn2vy = 7.957747154594767e-2*deltay*n/R;
	double 	dphi2_by_dn2vy = -1.0*dn1v_dot_n2v_by_dn2vy/(1 + -1.0*n3);
	double 	dn2vsqr_by_dn2vy = 2.0*n2vy;
	double 	dphi3_by_dn2vy = dn2vsqr_by_dn2vy*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphitot_by_dn2vy = dphi3_by_dn2vy + dphi2_by_dn2vy;
	double 	dn1v_dot_n2v_by_dn2vz = 7.957747154594767e-2*deltaz*n/R;
	double 	dphi2_by_dn2vz = -1.0*dn1v_dot_n2v_by_dn2vz/(1 + -1.0*n3);
	double 	dn2vsqr_by_dn2vz = 2.0*n2vz;
	double 	dphi3_by_dn2vz = dn2vsqr_by_dn2vz*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphitot_by_dn2vz = dphi3_by_dn2vz + dphi2_by_dn2vz;
	double 	n0 = 7.957747154594767e-2*deltak*n/(R*R);
	double 	dphi1_by_dn3 = n0/(1 + -1.0*n3);
	double 	n1vx = 7.957747154594767e-2*deltax*n/R;
	double 	n1vy = 7.957747154594767e-2*deltay*n/R;
	double 	n1vz = 7.957747154594767e-2*deltaz*n/R;
	double 	n1v_dot_n2v = n1vz*n2vz + n1vy*n2vy + n1vx*n2vx;
	double 	dphi2_by_dn3 = (n1*n2 + -1.0*n1v_dot_n2v)/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dphi3_by_dn3 = n2*(n2vsqr*(-5.305164769729845e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((-5.305164769729845e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(0.1061032953945969*1/(1 + -1.0*n3) + 5.305164769729845e-2))/(1 + -1.0*n3) + (2.6525823848649224e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(1 + -1.0*n3) + 5.305164769729845e-2) + 2.6525823848649224e-2)/n3)/n3) + n2*n2*(1.768388256576615e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((1.768388256576615e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(-3.53677651315323e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2))/(1 + -1.0*n3) + (-8.841941282883075e-3*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(1 + -1.0*n3) + -1.768388256576615e-2) + -8.841941282883075e-3)/n3)/n3))/(1 + -1.0*n3);
	double 	dphitot_by_dn3 = dphi3_by_dn3 + dphi2_by_dn3 + dphi1_by_dn3;
	double 	dphitot_by_dn0 = -1.0*log(1 + -1.0*n3);
	double 	dphitot_by_dn1 = n2/(1 + -1.0*n3);
	double 	dn1v_dot_n2v_by_dn1vx = n2vx;
	double 	dphitot_by_dn1vx = -1.0*dn1v_dot_n2v_by_dn1vx/(1 + -1.0*n3);
	double 	dn1v_dot_n2v_by_dn1vy = n2vy;
	double 	dphitot_by_dn1vy = -1.0*dn1v_dot_n2v_by_dn1vy/(1 + -1.0*n3);
	double 	dn1v_dot_n2v_by_dn1vz = n2vz;
	double 	dphitot_by_dn1vz = -1.0*dn1v_dot_n2v_by_dn1vz/(1 + -1.0*n3);
	double 	deltaprime = -25.132741228718345*R;
	double 	deltaprimex = 0;
	double 	deltaprimey = 0;
	double 	deltaprimez = 0;
	double 	dAdR = kT*n*(deltaprimez*dphitot_by_dn2vz + deltaprimey*dphitot_by_dn2vy + deltaprimex*dphitot_by_dn2vx + -1.0*deltaprime*dphitot_by_dn2 + deltak*dphitot_by_dn3 + (dphitot_by_dn1vz*(7.957747154594767e-2*deltaz/R + 7.957747154594767e-2*deltaprimez) + dphitot_by_dn1vy*(7.957747154594767e-2*deltay/R + 7.957747154594767e-2*deltaprimey) + dphitot_by_dn1vx*(7.957747154594767e-2*deltax/R + 7.957747154594767e-2*deltaprimex) + dphitot_by_dn1*(-7.957747154594767e-2*deltak/R + -7.957747154594767e-2*deltaprime) + dphitot_by_dn0*(-0.15915494309189535*deltak/R + -7.957747154594767e-2*deltaprime)/R)/R);
	double 	deltak2 = 50.26548245743669*R*R;
	double 	gSigmaA = dAdR/(deltak2*kT*n*n);
	double 	eta_d = 4.1887902047863905*R*R*R*n;
	double 	c1 = lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + 2.25855;
	double 	c2 = lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + -0.66927;
	double 	c3 = lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576;
	double 	eta_eff = eta_d*(eta_d*(c3*eta_d + c2) + c1);
	double 	ghs = (1 + -0.5*eta_eff)/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff));
	double 	deta_eff_by_deta_d = eta_d*(3.0*c3*eta_d + 2.0*c2) + c1;
	double 	dghs_by_deta_d = deta_eff_by_deta_d*((3.0*1 + -1.5*eta_eff)/(1 + -1.0*eta_eff) + -0.5)/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff));
	double 	da1_by_deta_d = epsilon_dispersion*(dghs_by_deta_d*eta_d + ghs)*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
	double 	dc1_by_dlambda_dispersion = -1.50349*1 + 0.498868*lambda_dispersion;
	double 	dc2_by_dlambda_dispersion = 1.40049*1 + -1.655478*lambda_dispersion;
	double 	dc3_by_dlambda_dispersion = -15.0427*1 + 10.61654*lambda_dispersion;
	double 	deta_eff_by_dlambda_dispersion = eta_d*(eta_d*(dc3_by_dlambda_dispersion*eta_d + dc2_by_dlambda_dispersion) + dc1_by_dlambda_dispersion);
	double 	dghs_by_dlambda_dispersion = deta_eff_by_dlambda_dispersion*((3.0*1 + -1.5*eta_eff)/(1 + -1.0*eta_eff) + -0.5)/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff));
	double 	da1_by_dlambda_dispersion = epsilon_dispersion*eta_d*(-12.0*ghs*lambda_dispersion*lambda_dispersion + dghs_by_dlambda_dispersion*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0));
	double 	gSW = (-8.333333333333333e-2*da1_by_dlambda_dispersion*lambda_dispersion/eta_d + 0.25*da1_by_deta_d)/kT + gSigmaA;
	double 	deltasaft = boltz*gSW*kappa_association;
	double 	delta2k = 50.26548245743669*R*R;
	double 	nA = 1.9894367886486918e-2*delta2k*n/(R*R);
	double 	X = (0.25*sqrt(8.0*deltasaft*nA + 1) + -0.25)/(deltasaft*nA);
	double 	dn3_by_dx = step;
	double 	ddphitot_by_dn0_by_dx = dn3_by_dx/(1 + -1.0*n3);
	double 	dn2_by_dx = deltak;
	double 	ddphitot_by_dn1_by_dx = dn3_by_dx*n2/((1 + -1.0*n3)*(1 + -1.0*n3)) + dn2_by_dx/(1 + -1.0*n3);
	double 	ddn1v_dot_n2v_by_dn1vx_by_dx = deltax;
	double 	ddphitot_by_dn1vx_by_dx = -1.0*dn1v_dot_n2v_by_dn1vx*dn3_by_dx/((1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*ddn1v_dot_n2v_by_dn1vx_by_dx/(1 + -1.0*n3);
	double 	ddn1v_dot_n2v_by_dn1vy_by_dx = deltay;
	double 	ddphitot_by_dn1vy_by_dx = -1.0*dn1v_dot_n2v_by_dn1vy*dn3_by_dx/((1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*ddn1v_dot_n2v_by_dn1vy_by_dx/(1 + -1.0*n3);
	double 	ddn1v_dot_n2v_by_dn1vz_by_dx = deltaz;
	double 	ddphitot_by_dn1vz_by_dx = -1.0*dn1v_dot_n2v_by_dn1vz*dn3_by_dx/((1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*ddn1v_dot_n2v_by_dn1vz_by_dx/(1 + -1.0*n3);
	double 	dn1_by_dx = 7.957747154594767e-2*deltak/R;
	double 	ddphi2_by_dn2_by_dx = dn3_by_dx*n1/((1 + -1.0*n3)*(1 + -1.0*n3)) + dn1_by_dx/(1 + -1.0*n3);
	double 	dn2vx_by_dx = deltax;
	double 	dn2vy_by_dx = deltay;
	double 	dn2vz_by_dx = deltaz;
	double 	dn2vsqr_by_dx = 2.0*dn2vz_by_dx*n2vz + 2.0*dn2vy_by_dx*n2vy + 2.0*dn2vx_by_dx*n2vx;
	double 	ddphi3_by_dn2_by_dx = 2.0*dn3_by_dx*(n2vsqr*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3)) + n2*n2*((log(1 + -1.0*n3)*(2.6525823848649224e-2*1/n3 + -5.305164769729845e-2) + 2.6525823848649224e-2)/n3 + 2.6525823848649224e-2*log(1 + -1.0*n3)))/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn3_by_dx*n2vsqr/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2vsqr*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2vsqr*(log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn3_by_dx*n2vsqr*log(1 + -1.0*n3)/(n3*n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + -2.6525823848649224e-2*dn3_by_dx*n2*n2/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2*(2.6525823848649224e-2*1/n3 + -5.305164769729845e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2*(log(1 + -1.0*n3)*(2.6525823848649224e-2*1/n3 + -5.305164769729845e-2) + 2.6525823848649224e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + -2.6525823848649224e-2*dn3_by_dx*n2*n2*log(1 + -1.0*n3)/(n3*n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn2vsqr_by_dx*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3)) + 2.0*dn2_by_dx*n2*((log(1 + -1.0*n3)*(2.6525823848649224e-2*1/n3 + -5.305164769729845e-2) + 2.6525823848649224e-2)/n3 + 2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	ddphitot_by_dn2_by_dx = ddphi3_by_dn2_by_dx + ddphi2_by_dn2_by_dx;
	double 	ddn1v_dot_n2v_by_dn2vx_by_dx = 7.957747154594767e-2*deltax/R;
	double 	ddphi2_by_dn2vx_by_dx = -1.0*dn1v_dot_n2v_by_dn2vx*dn3_by_dx/((1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*ddn1v_dot_n2v_by_dn2vx_by_dx/(1 + -1.0*n3);
	double 	ddn2vsqr_by_dn2vx_by_dx = 2.0*dn2vx_by_dx;
	double 	ddphi3_by_dn2vx_by_dx = 2.0*dn2vsqr_by_dn2vx*dn3_by_dx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn2vsqr_by_dn2vx*dn3_by_dx*n2/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn2vsqr_by_dn2vx*dn3_by_dx*n2*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn2vsqr_by_dn2vx*dn3_by_dx*n2*(log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn2vsqr_by_dn2vx*dn3_by_dx*n2*log(1 + -1.0*n3)/(n3*n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn2_by_dx*dn2vsqr_by_dn2vx*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3)) + ddn2vsqr_by_dn2vx_by_dx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	ddphitot_by_dn2vx_by_dx = ddphi3_by_dn2vx_by_dx + ddphi2_by_dn2vx_by_dx;
	double 	ddn1v_dot_n2v_by_dn2vy_by_dx = 7.957747154594767e-2*deltay/R;
	double 	ddphi2_by_dn2vy_by_dx = -1.0*dn1v_dot_n2v_by_dn2vy*dn3_by_dx/((1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*ddn1v_dot_n2v_by_dn2vy_by_dx/(1 + -1.0*n3);
	double 	ddn2vsqr_by_dn2vy_by_dx = 2.0*dn2vy_by_dx;
	double 	ddphi3_by_dn2vy_by_dx = 2.0*dn2vsqr_by_dn2vy*dn3_by_dx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn2vsqr_by_dn2vy*dn3_by_dx*n2/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn2vsqr_by_dn2vy*dn3_by_dx*n2*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn2vsqr_by_dn2vy*dn3_by_dx*n2*(log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn2vsqr_by_dn2vy*dn3_by_dx*n2*log(1 + -1.0*n3)/(n3*n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn2_by_dx*dn2vsqr_by_dn2vy*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3)) + ddn2vsqr_by_dn2vy_by_dx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	ddphitot_by_dn2vy_by_dx = ddphi3_by_dn2vy_by_dx + ddphi2_by_dn2vy_by_dx;
	double 	ddn1v_dot_n2v_by_dn2vz_by_dx = 7.957747154594767e-2*deltaz/R;
	double 	ddphi2_by_dn2vz_by_dx = -1.0*dn1v_dot_n2v_by_dn2vz*dn3_by_dx/((1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*ddn1v_dot_n2v_by_dn2vz_by_dx/(1 + -1.0*n3);
	double 	ddn2vsqr_by_dn2vz_by_dx = 2.0*dn2vz_by_dx;
	double 	ddphi3_by_dn2vz_by_dx = 2.0*dn2vsqr_by_dn2vz*dn3_by_dx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn2vsqr_by_dn2vz*dn3_by_dx*n2/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn2vsqr_by_dn2vz*dn3_by_dx*n2*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn2vsqr_by_dn2vz*dn3_by_dx*n2*(log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn2vsqr_by_dn2vz*dn3_by_dx*n2*log(1 + -1.0*n3)/(n3*n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn2_by_dx*dn2vsqr_by_dn2vz*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3)) + ddn2vsqr_by_dn2vz_by_dx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	ddphitot_by_dn2vz_by_dx = ddphi3_by_dn2vz_by_dx + ddphi2_by_dn2vz_by_dx;
	double 	dn0_by_dx = 7.957747154594767e-2*deltak/(R*R);
	double 	ddphi1_by_dn3_by_dx = dn3_by_dx*n0/((1 + -1.0*n3)*(1 + -1.0*n3)) + dn0_by_dx/(1 + -1.0*n3);
	double 	dn1vx_by_dx = 7.957747154594767e-2*deltax/R;
	double 	dn1vy_by_dx = 7.957747154594767e-2*deltay/R;
	double 	dn1vz_by_dx = 7.957747154594767e-2*deltaz/R;
	double 	dn1v_dot_n2v_by_dx = dn2vz_by_dx*n1vz + dn2vy_by_dx*n1vy + dn2vx_by_dx*n1vx + dn1vz_by_dx*n2vz + dn1vy_by_dx*n2vy + dn1vx_by_dx*n2vx;
	double 	ddphi2_by_dn3_by_dx = 2.0*dn3_by_dx*(n1*n2 + -1.0*n1v_dot_n2v)/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn2_by_dx*n1/((1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn1v_dot_n2v_by_dx/((1 + -1.0*n3)*(1 + -1.0*n3)) + dn1_by_dx*n2/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	ddphi3_by_dn3_by_dx = -1.768388256576615e-2*dn3_by_dx*n2*n2*n2/(((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + 3.53677651315323e-2*dn3_by_dx*n2*n2*n2*log(1 + -1.0*n3)/(((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + dn3_by_dx*n2*n2*n2*(1.768388256576615e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(-3.53677651315323e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2))/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2*n2*(-3.53677651315323e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + 1.768388256576615e-2*dn3_by_dx*n2*n2*n2/(n3*((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + -3.53677651315323e-2*dn3_by_dx*n2*n2*n2*log(1 + -1.0*n3)/(n3*((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + -1.0*dn3_by_dx*n2*n2*n2*((1.768388256576615e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(-3.53677651315323e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2))/(1 + -1.0*n3) + (-8.841941282883075e-3*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(1 + -1.0*n3) + -1.768388256576615e-2) + -8.841941282883075e-3)/n3)/(n3*n3*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2*n2*((1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(1 + -1.0*n3) + -1.768388256576615e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + -8.841941282883075e-3*dn3_by_dx*n2*n2*n2/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn3_by_dx*n2*n2*n2*log(1 + -1.0*n3)*(1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + 1.768388256576615e-2*dn3_by_dx*n2*n2*n2*log(1 + -1.0*n3)/(n3*n3*((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + -1.0*dn3_by_dx*n2*n2*n2*(-8.841941282883075e-3*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(1 + -1.0*n3) + -1.768388256576615e-2) + -8.841941282883075e-3)/(n3*n3*n3*(1 + -1.0*n3)) + 1.768388256576615e-2*dn3_by_dx*n2*n2*n2*log(1 + -1.0*n3)/((n3*n3)*(n3*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn3_by_dx*n2*(n2vsqr*(-5.305164769729845e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((-5.305164769729845e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(0.1061032953945969*1/(1 + -1.0*n3) + 5.305164769729845e-2))/(1 + -1.0*n3) + (2.6525823848649224e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(1 + -1.0*n3) + 5.305164769729845e-2) + 2.6525823848649224e-2)/n3)/n3) + n2*n2*(1.768388256576615e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((1.768388256576615e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(-3.53677651315323e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2))/(1 + -1.0*n3) + (-8.841941282883075e-3*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(1 + -1.0*n3) + -1.768388256576615e-2) + -8.841941282883075e-3)/n3)/n3))/((1 + -1.0*n3)*(1 + -1.0*n3)) + 5.305164769729845e-2*dn3_by_dx*n2*n2vsqr/(((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + -0.1061032953945969*dn3_by_dx*n2*n2vsqr*log(1 + -1.0*n3)/(((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + dn3_by_dx*n2*n2vsqr*(-5.305164769729845e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(0.1061032953945969*1/(1 + -1.0*n3) + 5.305164769729845e-2))/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2vsqr*(0.1061032953945969*1/(1 + -1.0*n3) + 5.305164769729845e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -5.305164769729845e-2*dn3_by_dx*n2*n2vsqr/(n3*((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + 0.1061032953945969*dn3_by_dx*n2*n2vsqr*log(1 + -1.0*n3)/(n3*((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + -1.0*dn3_by_dx*n2*n2vsqr*((-5.305164769729845e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(0.1061032953945969*1/(1 + -1.0*n3) + 5.305164769729845e-2))/(1 + -1.0*n3) + (2.6525823848649224e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(1 + -1.0*n3) + 5.305164769729845e-2) + 2.6525823848649224e-2)/n3)/(n3*n3*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2vsqr*((-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(1 + -1.0*n3) + 5.305164769729845e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn3_by_dx*n2*n2vsqr/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn3_by_dx*n2*n2vsqr*log(1 + -1.0*n3)*(-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -5.305164769729845e-2*dn3_by_dx*n2*n2vsqr*log(1 + -1.0*n3)/(n3*n3*((1 + -1.0*n3)*(1 + -1.0*n3))*((1 + -1.0*n3)*(1 + -1.0*n3))) + -1.0*dn3_by_dx*n2*n2vsqr*(2.6525823848649224e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(1 + -1.0*n3) + 5.305164769729845e-2) + 2.6525823848649224e-2)/(n3*n3*n3*(1 + -1.0*n3)) + -5.305164769729845e-2*dn3_by_dx*n2*n2vsqr*log(1 + -1.0*n3)/((n3*n3)*(n3*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn2vsqr_by_dx*n2*(-5.305164769729845e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((-5.305164769729845e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(0.1061032953945969*1/(1 + -1.0*n3) + 5.305164769729845e-2))/(1 + -1.0*n3) + (2.6525823848649224e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(1 + -1.0*n3) + 5.305164769729845e-2) + 2.6525823848649224e-2)/n3)/n3)/(1 + -1.0*n3) + dn2_by_dx*(n2vsqr*(-5.305164769729845e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((-5.305164769729845e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(0.1061032953945969*1/(1 + -1.0*n3) + 5.305164769729845e-2))/(1 + -1.0*n3) + (2.6525823848649224e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((-5.305164769729845e-2*1/(1 + -1.0*n3) + 5.305164769729845e-2*1/n3 + -0.1061032953945969)/(1 + -1.0*n3) + 5.305164769729845e-2) + 2.6525823848649224e-2)/n3)/n3) + n2*n2*(1.768388256576615e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((1.768388256576615e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(-3.53677651315323e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2))/(1 + -1.0*n3) + (-8.841941282883075e-3*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(1 + -1.0*n3) + -1.768388256576615e-2) + -8.841941282883075e-3)/n3)/n3))/(1 + -1.0*n3) + 2.0*dn2_by_dx*n2*n2*(1.768388256576615e-2*log(1 + -1.0*n3)/((1 + -1.0*n3)*(1 + -1.0*n3)) + ((1.768388256576615e-2*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*(-3.53677651315323e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2))/(1 + -1.0*n3) + (-8.841941282883075e-3*1/(1 + -1.0*n3) + log(1 + -1.0*n3)*((1.768388256576615e-2*1/(1 + -1.0*n3) + -1.768388256576615e-2*1/n3 + 3.53677651315323e-2)/(1 + -1.0*n3) + -1.768388256576615e-2) + -8.841941282883075e-3)/n3)/n3)/(1 + -1.0*n3);
	double 	ddphitot_by_dn3_by_dx = ddphi3_by_dn3_by_dx + ddphi2_by_dn3_by_dx + ddphi1_by_dn3_by_dx;
	double 	ddAdR_by_dx = kT*(deltaprimez*dphitot_by_dn2vz + deltaprimey*dphitot_by_dn2vy + deltaprimex*dphitot_by_dn2vx + -1.0*deltaprime*dphitot_by_dn2 + deltak*dphitot_by_dn3 + (dphitot_by_dn1vz*(7.957747154594767e-2*deltaz/R + 7.957747154594767e-2*deltaprimez) + dphitot_by_dn1vy*(7.957747154594767e-2*deltay/R + 7.957747154594767e-2*deltaprimey) + dphitot_by_dn1vx*(7.957747154594767e-2*deltax/R + 7.957747154594767e-2*deltaprimex) + dphitot_by_dn1*(-7.957747154594767e-2*deltak/R + -7.957747154594767e-2*deltaprime) + dphitot_by_dn0*(-0.15915494309189535*deltak/R + -7.957747154594767e-2*deltaprime)/R)/R) + ddphitot_by_dn3_by_dx*deltak*kT*n + ddphitot_by_dn2vz_by_dx*deltaprimez*kT*n + ddphitot_by_dn2vy_by_dx*deltaprimey*kT*n + ddphitot_by_dn2vx_by_dx*deltaprimex*kT*n + -1.0*ddphitot_by_dn2_by_dx*deltaprime*kT*n + ddphitot_by_dn1vz_by_dx*kT*n*(7.957747154594767e-2*deltaz/R + 7.957747154594767e-2*deltaprimez)/R + ddphitot_by_dn1vy_by_dx*kT*n*(7.957747154594767e-2*deltay/R + 7.957747154594767e-2*deltaprimey)/R + ddphitot_by_dn1vx_by_dx*kT*n*(7.957747154594767e-2*deltax/R + 7.957747154594767e-2*deltaprimex)/R + ddphitot_by_dn1_by_dx*kT*n*(-7.957747154594767e-2*deltak/R + -7.957747154594767e-2*deltaprime)/R + ddphitot_by_dn0_by_dx*kT*n*(-0.15915494309189535*deltak/R + -7.957747154594767e-2*deltaprime)/(R*R);
	double 	dgSigmaA_by_dx = ddAdR_by_dx/(deltak2*kT*n*n) + -2.0*dAdR/(deltak2*kT*n*n*n);
	double 	deta_d_by_dx = 4.1887902047863905*R*R*R;
	double 	ddeta_eff_by_deta_d_by_dx = deta_d_by_dx*(3.0*c3*eta_d + 2.0*c2) + 3.0*c3*deta_d_by_dx*eta_d;
	double 	deta_eff_by_dx = deta_d_by_dx*(eta_d*(c3*eta_d + c2) + c1) + deta_d_by_dx*eta_d*(c3*eta_d + c2) + c3*deta_d_by_dx*eta_d*eta_d;
	double 	ddghs_by_deta_d_by_dx = 3.0*deta_eff_by_deta_d*deta_eff_by_dx*((3.0*1 + -1.5*eta_eff)/(1 + -1.0*eta_eff) + -0.5)/(((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))) + -1.5*deta_eff_by_deta_d*deta_eff_by_dx/(((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))) + deta_eff_by_deta_d*deta_eff_by_dx*(3.0*1 + -1.5*eta_eff)/((1 + -1.0*eta_eff)*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))) + ddeta_eff_by_deta_d_by_dx*((3.0*1 + -1.5*eta_eff)/(1 + -1.0*eta_eff) + -0.5)/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff));
	double 	dghs_by_dx = -0.5*deta_eff_by_dx/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)) + 3.0*deta_eff_by_dx*(1 + -0.5*eta_eff)/(((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)));
	double 	dda1_by_deta_d_by_dx = dghs_by_dx*epsilon_dispersion*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0) + deta_d_by_dx*dghs_by_deta_d*epsilon_dispersion*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0) + ddghs_by_deta_d_by_dx*epsilon_dispersion*eta_d*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
	double 	ddeta_eff_by_dlambda_dispersion_by_dx = deta_d_by_dx*(eta_d*(dc3_by_dlambda_dispersion*eta_d + dc2_by_dlambda_dispersion) + dc1_by_dlambda_dispersion) + deta_d_by_dx*eta_d*(dc3_by_dlambda_dispersion*eta_d + dc2_by_dlambda_dispersion) + dc3_by_dlambda_dispersion*deta_d_by_dx*eta_d*eta_d;
	double 	ddghs_by_dlambda_dispersion_by_dx = 3.0*deta_eff_by_dlambda_dispersion*deta_eff_by_dx*((3.0*1 + -1.5*eta_eff)/(1 + -1.0*eta_eff) + -0.5)/(((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))) + -1.5*deta_eff_by_dlambda_dispersion*deta_eff_by_dx/(((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))) + deta_eff_by_dlambda_dispersion*deta_eff_by_dx*(3.0*1 + -1.5*eta_eff)/((1 + -1.0*eta_eff)*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))*((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff))) + ddeta_eff_by_dlambda_dispersion_by_dx*((3.0*1 + -1.5*eta_eff)/(1 + -1.0*eta_eff) + -0.5)/((1 + -1.0*eta_eff)*(1 + -1.0*eta_eff)*(1 + -1.0*eta_eff));
	double 	dda1_by_dlambda_dispersion_by_dx = -12.0*dghs_by_dx*epsilon_dispersion*eta_d*lambda_dispersion*lambda_dispersion + deta_d_by_dx*epsilon_dispersion*(-12.0*ghs*lambda_dispersion*lambda_dispersion + dghs_by_dlambda_dispersion*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0)) + ddghs_by_dlambda_dispersion_by_dx*epsilon_dispersion*eta_d*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
	double 	dgSW_by_dx = -8.333333333333333e-2*dda1_by_dlambda_dispersion_by_dx*lambda_dispersion/(eta_d*kT) + 0.25*dda1_by_deta_d_by_dx/kT + 8.333333333333333e-2*da1_by_dlambda_dispersion*deta_d_by_dx*lambda_dispersion/(eta_d*eta_d*kT) + dgSigmaA_by_dx;
	double 	ddeltasaft_by_dx = boltz*dgSW_by_dx*kappa_association;
	double 	dnA_by_dx = 1.9894367886486918e-2*delta2k/(R*R);
	double 	dX_by_dx = dnA_by_dx/(nA*sqrt(8.0*deltasaft*nA + 1)) + -1.0*dnA_by_dx*(0.25*sqrt(8.0*deltasaft*nA + 1) + -0.25)/(deltasaft*nA*nA) + ddeltasaft_by_dx/(deltasaft*sqrt(8.0*deltasaft*nA + 1)) + -1.0*ddeltasaft_by_dx*(0.25*sqrt(8.0*deltasaft*nA + 1) + -0.25)/(deltasaft*deltasaft*nA);
	double 	dFassoc_by_dx = kT*(2.0*1 + 4.0*log(X) + -2.0*X) + -2.0*dX_by_dx*kT*n + 4.0*dX_by_dx*kT*n/X;
	double 	a1 = epsilon_dispersion*eta_d*ghs*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
	double 	da1_by_dx = dghs_by_dx*epsilon_dispersion*eta_d*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0) + deta_d_by_dx*epsilon_dispersion*ghs*(-4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + 4.0);
	double 	da1integrated_by_dx = da1_by_dx*n + a1;
	double 	KHS = (eta_d*(eta_d*(eta_d*(-4.0*1 + eta_d) + 6.0) + -4.0) + 1)/(eta_d*(4.0*1 + 4.0*eta_d) + 1);
	double 	a2 = 0.5*KHS*da1_by_deta_d*epsilon_dispersion*eta_d;
	double 	dKHS_by_dx = deta_d_by_dx*(eta_d*(eta_d*(-4.0*1 + eta_d) + 6.0) + -4.0)/(eta_d*(4.0*1 + 4.0*eta_d) + 1) + -1.0*deta_d_by_dx*(4.0*1 + 4.0*eta_d)*(eta_d*(eta_d*(eta_d*(-4.0*1 + eta_d) + 6.0) + -4.0) + 1)/((eta_d*(4.0*1 + 4.0*eta_d) + 1)*(eta_d*(4.0*1 + 4.0*eta_d) + 1)) + deta_d_by_dx*eta_d*eta_d*eta_d/(eta_d*(4.0*1 + 4.0*eta_d) + 1) + deta_d_by_dx*eta_d*eta_d*(-4.0*1 + eta_d)/(eta_d*(4.0*1 + 4.0*eta_d) + 1) + deta_d_by_dx*eta_d*(eta_d*(-4.0*1 + eta_d) + 6.0)/(eta_d*(4.0*1 + 4.0*eta_d) + 1) + -4.0*deta_d_by_dx*eta_d*(eta_d*(eta_d*(eta_d*(-4.0*1 + eta_d) + 6.0) + -4.0) + 1)/((eta_d*(4.0*1 + 4.0*eta_d) + 1)*(eta_d*(4.0*1 + 4.0*eta_d) + 1));
	double 	da2_by_dx = 0.5*dKHS_by_dx*da1_by_deta_d*epsilon_dispersion*eta_d + 0.5*KHS*dda1_by_deta_d_by_dx*epsilon_dispersion*eta_d + 0.5*KHS*da1_by_deta_d*deta_d_by_dx*epsilon_dispersion;
	double 	da2integrated_by_dx = da2_by_dx*n/kT + a2/kT;
	double 	dFdisp_by_dx = da2integrated_by_dx + da1integrated_by_dx;
	double 	gpermol = 1822.8885;
	double 	mH2O = 18.01528*gpermol;
	double 	dFideal_by_dx = kT*(-1.0*1 + log(15.74960994572242*n/(sqrt(kT)*kT*sqrt(mH2O)*mH2O))) + kT;
	double 	dphi1_by_dx = dn3_by_dx*n0/(1 + -1.0*n3) + -1.0*dn0_by_dx*log(1 + -1.0*n3);
	double 	dphi2_by_dx = dn3_by_dx*(n1*n2 + -1.0*n1v_dot_n2v)/((1 + -1.0*n3)*(1 + -1.0*n3)) + dn2_by_dx*n1/(1 + -1.0*n3) + -1.0*dn1v_dot_n2v_by_dx/(1 + -1.0*n3) + dn1_by_dx*n2/(1 + -1.0*n3);
	double 	dphi3_by_dx = -8.841941282883075e-3*dn3_by_dx*n2*n2*n2/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2*n2*(8.841941282883075e-3*1/n3 + -1.768388256576615e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2*n2*(log(1 + -1.0*n3)*(8.841941282883075e-3*1/n3 + -1.768388256576615e-2) + 8.841941282883075e-3)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + -8.841941282883075e-3*dn3_by_dx*n2*n2*n2*log(1 + -1.0*n3)/(n3*n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.0*dn3_by_dx*n2*(n2vsqr*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3)) + n2*n2*((log(1 + -1.0*n3)*(8.841941282883075e-3*1/n3 + -1.768388256576615e-2) + 8.841941282883075e-3)/n3 + 8.841941282883075e-3*log(1 + -1.0*n3)))/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn3_by_dx*n2*n2vsqr/((1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2vsqr*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2)/(n3*(1 + -1.0*n3)*(1 + -1.0*n3)*(1 + -1.0*n3)) + -1.0*dn3_by_dx*n2*n2vsqr*(log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/(n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + 2.6525823848649224e-2*dn3_by_dx*n2*n2vsqr*log(1 + -1.0*n3)/(n3*n3*n3*(1 + -1.0*n3)*(1 + -1.0*n3)) + dn2vsqr_by_dx*n2*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3)) + dn2_by_dx*(n2vsqr*((log(1 + -1.0*n3)*(-2.6525823848649224e-2*1/n3 + 5.305164769729845e-2) + -2.6525823848649224e-2)/n3 + -2.6525823848649224e-2*log(1 + -1.0*n3)) + n2*n2*((log(1 + -1.0*n3)*(8.841941282883075e-3*1/n3 + -1.768388256576615e-2) + 8.841941282883075e-3)/n3 + 8.841941282883075e-3*log(1 + -1.0*n3)))/((1 + -1.0*n3)*(1 + -1.0*n3)) + 2.0*dn2_by_dx*n2*n2*((log(1 + -1.0*n3)*(8.841941282883075e-3*1/n3 + -1.768388256576615e-2) + 8.841941282883075e-3)/n3 + 8.841941282883075e-3*log(1 + -1.0*n3))/((1 + -1.0*n3)*(1 + -1.0*n3));
	double 	dwhitebear_by_dx = dphi3_by_dx*kT + dphi2_by_dx*kT + dphi1_by_dx*kT;
	double 	dFSAFT_by_dx = mu + dwhitebear_by_dx + dFideal_by_dx + dFdisp_by_dx + dFassoc_by_dx;
	output = dFSAFT_by_dx;

	return output;

}

double d_by_dT(double , double ) const {
	assert(0); // fail
	return 0;

}

Functional grad(const Functional &ingrad, const Functional &x, bool ) const {
	return ingrad;}

Functional grad_T(const Functional &ingradT) const {
	return ingradT;}

void grad(const GridDescription &gd, double kT, const VectorXd &x, const VectorXd &ingrad, VectorXd *outgrad, VectorXd *outpgrad) const {
	if (oldkT != kT) {
		oldkT = kT;

	}
	VectorXcd ktemp0(gd.NxNyNzOver2);
	ktemp0 = fft(gd, x);

	VectorXcd ktemp1(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp1[0] = ktemp0[i]*(12.566370614359172*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp1[i] = ktemp0[i]*(12.566370614359172*R*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1);
	}

	VectorXd rtemp2(gd.NxNyNz);
	rtemp2 = ifft(gd, ktemp1);

	VectorXcd ktemp3(gd.NxNyNzOver2);
	ktemp3[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp3[i] = ktemp0[i]*(12.566370614359172*complex(0,1)*k_i[0]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1);
	}

	VectorXd rtemp4(gd.NxNyNz);
	rtemp4 = ifft(gd, ktemp3);

	VectorXcd ktemp5(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp5[0] = ktemp0[i]*(4.188790204786391*R*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp5[i] = ktemp0[i]*(12.566370614359172*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*t2*cos(t2) + sin(t2))/(sqrt(t1)*t1));
	}

	VectorXd rtemp6(gd.NxNyNz);
	rtemp6 = ifft(gd, ktemp5);

	ktemp5.resize(0); // KSpace
	VectorXd rtemp7(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp7[i] = rtemp2[i]/(1 + -1.0*rtemp6[i]);
	}

	VectorXcd ktemp8(gd.NxNyNzOver2);
	ktemp8 = fft(gd, rtemp7);

	VectorXcd ktemp9(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp9[0] = ktemp0[i]*(50.26548245743669*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp9[i] = ktemp0[i]*(25.132741228718345*R*sin(2.0*R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1);
	}

	VectorXd rtemp10(gd.NxNyNz);
	rtemp10 = ifft(gd, ktemp9);

	ktemp9.resize(0); // KSpace
	VectorXcd ktemp11(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp11[0] = ktemp0[i]*(50.26548245743669*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp11[i] = ktemp0[i]*(25.132741228718345*R*sin(2.0*R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1);
	}

	VectorXd rtemp12(gd.NxNyNz);
	rtemp12 = ifft(gd, ktemp11);

	ktemp11.resize(0); // KSpace
	VectorXcd ktemp13(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp13[0] = ktemp0[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp13[i] = ktemp0[i]*exp(-4.0*R*R*lambda_dispersion*lambda_dispersion*length_scaling*length_scaling*(0.5*k_i[2]*k_i[2] + 0.5*k_i[1]*k_i[1] + 0.5*k_i[0]*k_i[0]));
	}

	VectorXd rtemp14(gd.NxNyNz);
	rtemp14 = ifft(gd, ktemp13);

	ktemp13.resize(0); // KSpace
	VectorXd rtemp15(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp15[i] = -1.0*rtemp4[i]/(1 + -1.0*rtemp6[i]);
	}

	VectorXcd ktemp16(gd.NxNyNzOver2);
	ktemp16 = fft(gd, rtemp15);

	rtemp15.resize(0); // Realspace
	VectorXcd ktemp17(gd.NxNyNzOver2);
	ktemp17[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp17[i] = ktemp0[i]*(12.566370614359172*complex(0,1)*k_i[1]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1);
	}

	VectorXd rtemp18(gd.NxNyNz);
	rtemp18 = ifft(gd, ktemp17);

	VectorXd rtemp19(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp19[i] = -1.0*rtemp18[i]/(1 + -1.0*rtemp6[i]);
	}

	VectorXcd ktemp20(gd.NxNyNzOver2);
	ktemp20 = fft(gd, rtemp19);

	rtemp19.resize(0); // Realspace
	ktemp0[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp0[i] = ktemp0[i]*(12.566370614359172*complex(0,1)*k_i[2]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1);
	}

	VectorXd rtemp22(gd.NxNyNz);
	rtemp22 = ifft(gd, ktemp0);

	VectorXd rtemp23(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp23[i] = -1.0*rtemp22[i]/(1 + -1.0*rtemp6[i]);
	}

	VectorXcd ktemp24(gd.NxNyNzOver2);
	ktemp24 = fft(gd, rtemp23);

	rtemp23.resize(0); // Realspace
	VectorXd rtemp25(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp25[i] = -1.0*log(1 + -1.0*rtemp6[i]);
	}

	VectorXcd ktemp26(gd.NxNyNzOver2);
	ktemp26 = fft(gd, rtemp25);

	VectorXd rtemp27(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		const double t2 = rtemp18[i]*rtemp18[i];
		const double t3 = rtemp22[i]*rtemp22[i];
		const double t4 = rtemp2[i]*rtemp2[i];
		const double t5 = rtemp4[i]*rtemp4[i];
		const double t6 = 1/t1;
		const double t7 = 1/rtemp6[i];
		rtemp27[i] = t7*t7*rtemp7[i]*(-3.0*t5 + t4 + -3.0*t3 + -3.0*t2)*(t6*((-1.768388256576615e-2*t7 + 1.768388256576615e-2*t6)*(log(t1)/(t6*t6) + rtemp6[i]) + 8.841941282883075e-3) + -8.841941282883075e-3*1 + 1.768388256576615e-2*rtemp25[i]) + t6*(t6*(-7.957747154594767e-2*t5 + 7.957747154594767e-2*t4 + -7.957747154594767e-2*t3 + -7.957747154594767e-2*t2) + 7.957747154594767e-2*rtemp2[i]/R)/R;
	}

	VectorXcd ktemp28(gd.NxNyNzOver2);
	ktemp28 = fft(gd, rtemp27);

	rtemp27.resize(0); // Realspace
	VectorXd rtemp29(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		rtemp29[i] = ((t1*t1*log(t1) + rtemp6[i])*(-2.6525823848649224e-2*rtemp4[i]*rtemp4[i] + 2.6525823848649224e-2*rtemp2[i]*rtemp2[i] + -2.6525823848649224e-2*rtemp22[i]*rtemp22[i] + -2.6525823848649224e-2*rtemp18[i]*rtemp18[i])/(t1*rtemp6[i]*rtemp6[i]) + 7.957747154594767e-2*rtemp2[i]/R)/t1;
	}

	VectorXcd ktemp30(gd.NxNyNzOver2);
	ktemp30 = fft(gd, rtemp29);

	rtemp29.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp30[0] = -1.0*ktemp30[i]*(-25.132741228718345*R) + ktemp28[i]*(12.566370614359172*R*R);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sin(R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t2 = exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t3 = 1/(sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t4 = t1*t3;
		ktemp30[i] = -1.0*ktemp30[i]*(12.566370614359172*t2*(-1.0*R*cos(R/t3) + -1.0*t4)) + ktemp28[i]*(12.566370614359172*R*t2*t4);
	}

	ktemp28.resize(0); // KSpace
	VectorXd rtemp32(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		rtemp32[i] = rtemp4[i]*(-5.305164769729845e-2*rtemp7[i]*(t1*t1*log(t1) + rtemp6[i])/(rtemp6[i]*rtemp6[i]) + -7.957747154594767e-2*1/R)/t1;
	}

	VectorXcd ktemp33(gd.NxNyNzOver2);
	ktemp33 = fft(gd, rtemp32);

	rtemp32.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp33[0] = ktemp30[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp33[i] = ktemp33[i]*(12.566370614359172*R*complex(0,1)*k_i[0]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1) + ktemp30[i];
	}

	ktemp30.resize(0); // KSpace
	VectorXd rtemp35(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		rtemp35[i] = rtemp18[i]*(-5.305164769729845e-2*rtemp7[i]*(t1*t1*log(t1) + rtemp6[i])/(rtemp6[i]*rtemp6[i]) + -7.957747154594767e-2*1/R)/t1;
	}

	VectorXcd ktemp36(gd.NxNyNzOver2);
	ktemp36 = fft(gd, rtemp35);

	rtemp35.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp36[0] = ktemp33[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp36[i] = ktemp36[i]*(12.566370614359172*R*complex(0,1)*k_i[1]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1) + ktemp33[i];
	}

	ktemp33.resize(0); // KSpace
	VectorXd rtemp38(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		rtemp38[i] = rtemp22[i]*(-5.305164769729845e-2*rtemp7[i]*(t1*t1*log(t1) + rtemp6[i])/(rtemp6[i]*rtemp6[i]) + -7.957747154594767e-2*1/R)/t1;
	}

	VectorXcd ktemp39(gd.NxNyNzOver2);
	ktemp39 = fft(gd, rtemp38);

	rtemp38.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp39[0] = ktemp36[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp39[i] = ktemp39[i]*(12.566370614359172*R*complex(0,1)*k_i[2]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1) + ktemp36[i];
	}

	ktemp36.resize(0); // KSpace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp39[0] = (-7.957747154594767e-2*ktemp8[i]*(-25.132741228718345*R) + (-7.957747154594767e-2*ktemp8[i]*(12.566370614359172*R*R) + ktemp26[i]*(-0.15915494309189535*(12.566370614359172*R*R)/R + -7.957747154594767e-2*(-25.132741228718345*R)))/R)/R + ktemp39[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t2 = -1.0*sin(R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]))/(sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0])) + R*cos(R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t3 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t4 = 12.566370614359172*t1*(-1.0*sin(R*sqrt(t3))/(sqrt(t3)) + -1.0*R*cos(R*sqrt(t3)));
		const complex t5 = 12.566370614359172*R*t1*sin(R*sqrt(t3))/(sqrt(t3));
		const complex t6 = sin(R*sqrt(t3));
		ktemp39[i] = (7.957747154594767e-2*ktemp24[i]*(12.566370614359172*R*t1*t6*complex(0,1)*k_i[2]/(sqrt(t3))) + 7.957747154594767e-2*ktemp20[i]*(12.566370614359172*R*t1*t6*complex(0,1)*k_i[1]/(sqrt(t3))) + 7.957747154594767e-2*ktemp16[i]*(12.566370614359172*R*t1*t6*complex(0,1)*k_i[0]/(sqrt(t3))) + -7.957747154594767e-2*t4*ktemp8[i] + (-7.957747154594767e-2*ktemp26[i]*(2.0*t5/R + t4) + 7.957747154594767e-2*ktemp24[i]*(12.566370614359172*t1*t2*complex(0,1)*k_i[2]/t3) + 7.957747154594767e-2*ktemp20[i]*(12.566370614359172*t1*t2*complex(0,1)*k_i[1]/t3) + 7.957747154594767e-2*ktemp16[i]*(12.566370614359172*t1*t2*complex(0,1)*k_i[0]/t3) + -7.957747154594767e-2*t5*ktemp8[i])/R)/R + ktemp39[i];
	}

	ktemp26.resize(0); // KSpace
	ktemp24.resize(0); // KSpace
	ktemp20.resize(0); // KSpace
	ktemp16.resize(0); // KSpace
	ktemp8.resize(0); // KSpace
	VectorXd rtemp42(gd.NxNyNz);
	rtemp42 = ifft(gd, ktemp39);

	ktemp39.resize(0); // KSpace
	VectorXd rtemp43(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp43[i] = R*kT*t4*x[i]*(23812.820490470258*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -3789.9280900183135*1 + -1894.9640450091567*t3)/rtemp12[i] + 301.59289474462014*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp44(gd.NxNyNzOver2);
	ktemp44 = fft(gd, rtemp43);

	rtemp43.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp44[0] = R*ktemp44[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp44[i] = ktemp44[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp46(gd.NxNyNz);
	rtemp46 = ifft(gd, ktemp44);

	ktemp44.resize(0); // KSpace
	VectorXd rtemp47(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp47[i] = R*kT*t4*x[i]*(-210.55156055657295*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + 33.510321638291124*1 + 16.755160819145562*t3)/rtemp12[i] + -2.6666666666666665*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp48(gd.NxNyNzOver2);
	ktemp48 = fft(gd, rtemp47);

	rtemp47.resize(0); // Realspace
	VectorXcd ktemp49(gd.NxNyNzOver2);
	ktemp49[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp49[i] = ktemp48[i]*complex(0,1)*k_i[0]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp50(gd.NxNyNz);
	rtemp50 = ifft(gd, ktemp49);

	ktemp49.resize(0); // KSpace
	VectorXd rtemp51(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT;
		rtemp51[i] = kT*x[i]/(t1*rtemp12[i]*sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(-1.0*1 + exp(epsilon_association/kT))/(R*R) + 1));
	}

	VectorXcd ktemp52(gd.NxNyNzOver2);
	ktemp52 = fft(gd, rtemp51);

		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp52[0] = -2.0*R*ktemp52[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp52[i] = ktemp52[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp54(gd.NxNyNz);
	rtemp54 = ifft(gd, ktemp52);

	ktemp52.resize(0); // KSpace
	VectorXd rtemp55(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1);
		rtemp55[i] = 7.957747154594767e-2*kT*kappa_association*t1*rtemp10[i]*x[i]/(R*R*t2*rtemp12[i]*(-1.0*1 + t2));
	}

	VectorXcd ktemp56(gd.NxNyNzOver2);
	ktemp56 = fft(gd, rtemp55);

	rtemp55.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp56[0] = -2.0*R*ktemp56[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp56[i] = ktemp56[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp58(gd.NxNyNz);
	rtemp58 = ifft(gd, ktemp56);

	ktemp56.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp58[i] = -33.510321638291124*rtemp58[i] + 16.755160819145562*rtemp54[i];
	}

	rtemp54.resize(0); // Realspace
	VectorXd rtemp60(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1/(-1.0*1 + exp(epsilon_association/kT));
		const double t2 = rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT;
		rtemp60[i] = 50.26548245743669*R*R*kT*t1*x[i]*(sqrt(0.15915494309189535*kappa_association*t2*rtemp10[i]/(R*R*t1) + 1) + -1.0)/(kappa_association*t2*t2*rtemp10[i]*rtemp12[i]);
	}

	VectorXcd ktemp61(gd.NxNyNzOver2);
	ktemp61 = fft(gd, rtemp60);

	rtemp60.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp61[0] = -2.0*R*ktemp61[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp61[i] = ktemp61[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp63(gd.NxNyNz);
	rtemp63 = ifft(gd, ktemp61);

	ktemp61.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp63[i] = -4.1887902047863905*rtemp63[i] + rtemp58[i];
	}

	rtemp58.resize(0); // Realspace
	VectorXd rtemp65(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t2 = 4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0;
		const double t3 = R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion);
		const double t4 = 4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + t3;
		const double t5 = -2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t6 = R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776);
		const double t7 = lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion);
		const double t8 = 2.25855*1 + t7 + t6;
		rtemp65[i] = 4.0*kT*x[i]/(rtemp12[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*(t2*(-4.1887902047863905*R*R*R*rtemp14[i]*(3.0*t5*t8/t1 + -0.5*t8) + -1.0*t5)/(t1*t1*t1) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*t5 + -4.1887902047863905*R*R*R*t2*(3.0*t4*t5/t1 + -0.5*t4))/(t1*t1*t1))/kT));
	}

	VectorXcd ktemp66(gd.NxNyNzOver2);
	ktemp66 = fft(gd, rtemp65);

	rtemp65.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp66[0] = -2.0*R*ktemp66[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp66[i] = ktemp66[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp68(gd.NxNyNz);
	rtemp68 = ifft(gd, ktemp66);

	ktemp66.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp68[i] = 8.377580409572781*rtemp68[i] + rtemp63[i];
	}

	rtemp63.resize(0); // Realspace
	VectorXd rtemp70(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp70[i] = kT*x[i]*(t4*(631.6546816697189*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -100.53096491487338*1 + -50.26548245743669*t3) + 8.0*kappa_association*t1*t2*t3*rtemp10[i]/(R*R))/rtemp12[i];
	}

	VectorXcd ktemp71(gd.NxNyNzOver2);
	ktemp71 = fft(gd, rtemp70);

	rtemp70.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp71[0] = R*ktemp71[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp71[i] = ktemp71[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp73(gd.NxNyNz);
	rtemp73 = ifft(gd, ktemp71);

	VectorXd rtemp74(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp74[i] = kT*x[i]*(t4*(-631.6546816697189*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + 100.53096491487338*1 + 50.26548245743669*t3) + -8.0*kappa_association*t1*t2*t3*rtemp10[i]/(R*R))/rtemp12[i];
	}

	VectorXcd ktemp75(gd.NxNyNzOver2);
	ktemp75 = fft(gd, rtemp74);

	rtemp74.resize(0); // Realspace
	VectorXcd ktemp76(gd.NxNyNzOver2);
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp76[0] = 2.0*kT*ktemp1[i] + -2.0*R*ktemp75[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp76[i] = ktemp75[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2)) + 2.0*kT*ktemp1[i];
	}

	VectorXd rtemp77(gd.NxNyNz);
	rtemp77 = ifft(gd, ktemp76);

	ktemp76.resize(0); // KSpace
	VectorXd rtemp78(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp78[i] = -315.82734083485946*kT*t4*t4*x[i]/(kappa_association*t1*t2*rtemp10[i]*rtemp12[i]) + kT*x[i]*(t4*(50.26548245743669*1 + 25.132741228718345*t3) + -4.0*kappa_association*t1*t2*t3*rtemp10[i]/(R*R))/(R*R*rtemp12[i]);
	}

	VectorXcd ktemp79(gd.NxNyNzOver2);
	ktemp79 = fft(gd, rtemp78);

	rtemp78.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp3[0] = 2.0*kT*ktemp3[i]/R;
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		const complex t3 = sin(t2)/(sqrt(t1));
		ktemp3[i] = complex(0,1)*k_i[0]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(t3*ktemp75[i] + ktemp79[i]*(R*cos(t2) + -1.0*t3)/t1) + 2.0*kT*ktemp3[i]/R;
	}

	VectorXd rtemp81(gd.NxNyNz);
	rtemp81 = ifft(gd, ktemp3);

	ktemp3.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		const double t2 = t1*t1*log(t1) + rtemp6[i];
		const double t3 = 1/t1;
		const double t4 = 1/rtemp6[i];
		rtemp81[i] = rtemp4[i]*(t3*t3*rtemp73[i] + t2*t3*t3*t4*t4*rtemp68[i]) + rtemp2[i]*(t3*t4*t4*rtemp4[i]*(rtemp46[i]*(t3*(t2*(-1.768388256576615e-2*t4 + 1.768388256576615e-2*t3) + 8.841941282883075e-3) + -8.841941282883075e-3*1 + 1.768388256576615e-2*rtemp25[i]) + 0.6666666666666667*kT*t2*t3) + t2*t3*t3*t4*t4*rtemp50[i]) + t3*rtemp81[i];
	}

	VectorXcd ktemp83(gd.NxNyNzOver2);
	ktemp83 = fft(gd, rtemp81);

	rtemp81.resize(0); // Realspace
	VectorXcd ktemp84(gd.NxNyNzOver2);
	ktemp84[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp84[i] = ktemp48[i]*complex(0,1)*k_i[1]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp85(gd.NxNyNz);
	rtemp85 = ifft(gd, ktemp84);

	ktemp84.resize(0); // KSpace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp17[0] = 2.0*kT*ktemp17[i]/R;
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		const complex t3 = sin(t2)/(sqrt(t1));
		ktemp17[i] = complex(0,1)*k_i[1]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(t3*ktemp75[i] + ktemp79[i]*(R*cos(t2) + -1.0*t3)/t1) + 2.0*kT*ktemp17[i]/R;
	}

	VectorXd rtemp87(gd.NxNyNz);
	rtemp87 = ifft(gd, ktemp17);

	ktemp17.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		const double t2 = t1*t1*log(t1) + rtemp6[i];
		const double t3 = 1/t1;
		const double t4 = 1/rtemp6[i];
		rtemp87[i] = rtemp2[i]*(t3*t4*t4*rtemp18[i]*(rtemp46[i]*(t3*(t2*(-1.768388256576615e-2*t4 + 1.768388256576615e-2*t3) + 8.841941282883075e-3) + -8.841941282883075e-3*1 + 1.768388256576615e-2*rtemp25[i]) + 0.6666666666666667*kT*t2*t3) + t2*t3*t3*t4*t4*rtemp85[i]) + rtemp18[i]*(t3*t3*rtemp73[i] + t2*t3*t3*t4*t4*rtemp68[i]) + t3*rtemp87[i];
	}

	VectorXcd ktemp89(gd.NxNyNzOver2);
	ktemp89 = fft(gd, rtemp87);

	rtemp87.resize(0); // Realspace
	ktemp89[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp89[i] = ktemp89[i]*k_i[1] + ktemp83[i]*k_i[0];
	}

	ktemp83.resize(0); // KSpace
	ktemp48[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp48[i] = ktemp48[i]*complex(0,1)*k_i[2]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp92(gd.NxNyNz);
	rtemp92 = ifft(gd, ktemp48);

	ktemp48.resize(0); // KSpace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp75[0] = 2.0*kT*ktemp0[i]/R;
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		const complex t3 = sin(t2)/(sqrt(t1));
		ktemp75[i] = complex(0,1)*k_i[2]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(t3*ktemp75[i] + ktemp79[i]*(R*cos(t2) + -1.0*t3)/t1) + 2.0*kT*ktemp0[i]/R;
	}

	ktemp0.resize(0); // KSpace
	VectorXd rtemp94(gd.NxNyNz);
	rtemp94 = ifft(gd, ktemp75);

	ktemp75.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		const double t2 = t1*t1*log(t1) + rtemp6[i];
		const double t3 = 1/t1;
		const double t4 = 1/rtemp6[i];
		rtemp94[i] = rtemp2[i]*(t3*t4*t4*rtemp22[i]*(rtemp46[i]*(t3*(t2*(-1.768388256576615e-2*t4 + 1.768388256576615e-2*t3) + 8.841941282883075e-3) + -8.841941282883075e-3*1 + 1.768388256576615e-2*rtemp25[i]) + 0.6666666666666667*kT*t2*t3) + t2*t3*t3*t4*t4*rtemp92[i]) + rtemp22[i]*(t3*t3*rtemp73[i] + t2*t3*t3*t4*t4*rtemp68[i]) + t3*rtemp94[i];
	}

	rtemp73.resize(0); // Realspace
	rtemp46.resize(0); // Realspace
	VectorXcd ktemp96(gd.NxNyNzOver2);
	ktemp96 = fft(gd, rtemp94);

	rtemp94.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp96[0] = ktemp89[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp96[i] = ktemp96[i]*k_i[2] + ktemp89[i];
	}

	ktemp89.resize(0); // KSpace
	VectorXd rtemp98(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = -4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t4 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/t3 + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/(t3*t3*t3) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(-2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + 3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/t3 + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/(t3*t3*t3))/kT)/(R*R) + 1));
		const double t5 = -2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t6 = -50.26548245743669*lambda_dispersion*lambda_dispersion*t5 + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(-2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + 3.0*t5*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/t3 + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion));
		const double t7 = 4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0;
		const double t8 = -4.1887902047863905*R*R*R*rtemp14[i] + 1;
		const double t9 = 1/(4.0*(4.1887902047863905*R*R*R*rtemp14[i])*(4.1887902047863905*R*R*R*rtemp14[i]) + 16.755160819145562*R*R*R*rtemp14[i] + 1);
		const double t10 = lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855;
		const double t11 = lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion);
		const double t12 = t5*(10.1576*1 + t11)/t3;
		const double t13 = 10.1576*1 + t11;
		const double t14 = 3694.3299710673505*t5/t3 + -923.5824927668376;
		const double t15 = 4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion);
		const double t16 = rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*(t7*(-4.1887902047863905*R*R*R*rtemp14[i]*(3.0*t10*t5/t3 + -0.5*t10) + -1.0*t5)/(t3*t3*t3) + -7.957747154594767e-2*lambda_dispersion*t6/(t3*t3*t3))/kT;
		const double t17 = lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + -0.66927;
		const double t18 = R*R*R*rtemp14[i]*(52.63789013914324*R*R*R*t13*rtemp14[i] + 8.377580409572781*t17);
		const double t19 = lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion);
		const double t20 = t10*t5/t3;
		const double t21 = -587.9708763079075*1 + -293.98543815395374*t4;
		const double t22 = 11082.989913202051*t5/t3 + -1847.1649855336752;
		const double t23 = -15.0427*1 + 10.61654*lambda_dispersion;
		const double t24 = t13*t22/(t16*t16*t3*t3*t3);
		const double t25 = 46.78923567923842*t4*t5/t3 + -11.697308919809606*t4;
		const double t26 = R*R*R*rtemp14[i]*(17.54596337971441*R*R*R*t13*rtemp14[i] + 4.1887902047863905*t17);
		const double t27 = t5*(2.25855*1 + t26 + t19)/t3;
		const double t28 = t9*(-4.1887902047863905*R*R*R*rtemp14[i]*(3.0*t20 + -0.5*t10) + -1.0*t5);
		const double t29 = -0.5*t13*t21/t16 + 3.0000000000000004*t12*t21/t16;
		const double t30 = t15*t21/t16;
		const double t31 = t21/t16;
		const double t32 = t31*t5/t3;
		const double t33 = t10*t31;
		const double t34 = t29/(t3*t3*t3);
		const double t35 = 1322.9344716927922*t5/t3 + -220.48907861546533;
		const double t36 = 2.25855*1 + t26 + t19;
		const double t37 = 1.40049*1 + -1.655478*lambda_dispersion;
		const double t38 = 220.48907861546533*t5/t3 + -36.74817976924422;
		const double t39 = 140.3677070377153*t4*t5/t3 + -23.39461783961921*t4;
		const double t40 = t13*t39/(t3*t3*t3);
		const double t41 = -210.55156055657295*1 + -105.27578027828648*t4;
		const double t42 = 5.968310365946076e-2*t32 + -9.947183943243459e-3*t31;
		const double t43 = t41/t16;
		const double t44 = -1.50349*1 + 0.498868*lambda_dispersion;
		const double t45 = 3.0*t15*t5/t3 + -0.5*t15;
		const double t46 = t35*t36/(t16*t16);
		const double t47 = 3.0*t20 + -0.5*t10;
		const double t48 = -1.4248291449703753e-2*t31*t36 + 8.548974869822253e-2*t27*t31;
		const double t49 = lambda_dispersion*lambda_dispersion*t5;
		const double t50 = t36*t39;
		rtemp98[i] = epsilon_dispersion*lambda_dispersion*(1.3534253949046735e-4*t31*t6*x[i]/(t3*t3*t3*rtemp14[i]) + t49*x[i]*(-4.000000000000001*1 + -2.0*t4)/(t16*t3*t3*t3*rtemp14[i])) + R*(epsilon_dispersion*kappa_association*(t1*t2*t7*rtemp10[i]*x[i]*(-0.3333333333333333*t4*t47 + -7.124145724851877e-3*t50)/(t3*t3*t3) + lambda_dispersion*(rtemp10[i]*x[i]*(2.6525823848649224e-2*t1*t2*t4*t45*t7/(t3*t3*t3*rtemp14[i]) + 2.374715241617292e-3*t1*t2*t39*t44*t7/(t3*t3*t3)) + 2.8496582899407506e-2*lambda_dispersion*lambda_dispersion*t1*t2*t50*rtemp10[i]*x[i]/(t3*t3*t3))) + R*(epsilon_dispersion*lambda_dispersion*x[i]*(-2.094395102393195*R*R*R*t45*t7 + 7.105427357601002e-15*t49)/(kappa_association*t1*t16*t16*t2*t3*t3*t3*rtemp10[i]*rtemp14[i]) + R*(epsilon_dispersion*(t7*x[i]*((2.0943951023931953*epsilon_dispersion*t28*(t8*t8)*(t8*t8)/kT + -4.1887902047863905*t5 + -0.25*t48)/(t3*t3*t3) + -7.124145724851877e-3*t31*t47/(t3*t3*t3)) + lambda_dispersion*(x[i]*(0.11936620731892152*t42*t44*t7/(t3*t3*t3) + 5.66921503708585e-4*t31*t45*t7/(t3*t3*t3*rtemp14[i])) + lambda_dispersion*lambda_dispersion*t48*x[i]/(t3*t3*t3))) + R*(epsilon_dispersion*kappa_association*(t1*t2*t7*rtemp10[i]*rtemp14[i]*x[i]*(t17*(-12.566370614359172*t4*t5/t3 + 2.0943951023931953*t4)/(t3*t3*t3) + -0.3580986219567646*t10*t25*t36/((t3*t3)*(t3*t3))) + lambda_dispersion*(rtemp10[i]*x[i]*(1.9894367886486918e-2*t1*t2*t37*t39*t7*rtemp14[i]/(t3*t3*t3) + 2.8496582899407506e-2*t1*t15*t2*t25*t36*t7/((t3*t3)*(t3*t3))) + 0.11936620731892152*lambda_dispersion*lambda_dispersion*t1*t17*t2*t39*rtemp10[i]*rtemp14[i]*x[i]/(t3*t3*t3))) + R*(epsilon_dispersion*(t7*x[i]*(-26.31894506957162*t47/(t16*t16) + -5.968310365946076e-2*t46)/(t1*t2*t3*t3*t3*rtemp10[i]) + lambda_dispersion*(x[i]*(2.0943951023931953*t45*t7/(t1*t16*t16*t2*t3*t3*t3*rtemp14[i]) + 0.11936620731892152*t38*t44*t7/(t1*t16*t16*t2*t3*t3*t3))/rtemp10[i] + 0.23873241463784303*lambda_dispersion*lambda_dispersion*t46*x[i]/(t1*t2*t3*t3*t3*rtemp10[i])))/kappa_association + R*(epsilon_dispersion*(t7*rtemp14[i]*x[i]*(((-52.63789013914324*t36*t5 + -1.0*t20*t36*t43 + 0.25*t10*t36*t43)/t3 + 8.772981689857206*t36)/(t3*t3*t3) + t17*(0.125*t43 + -0.2685739664675735*t32)/(t3*t3*t3)) + lambda_dispersion*(x[i]*(t37*t42*t7*rtemp14[i]/(t3*t3*t3) + t7*(-7.124145724851877e-3*t30*t36 + 2.8496582899407506e-2*t27*t30)/((t3*t3)*(t3*t3))) + lambda_dispersion*lambda_dispersion*t17*rtemp14[i]*x[i]*(t43*t5/t3 + -5.968310365946076e-2*t31)/(t3*t3*t3)) + epsilon_dispersion*t7*rtemp14[i]*x[i]*((t8*t8)*(t8*t8)*t9*(19.81421779562699*1 + -26.31894506957162*t27 + 4.386490844928603*t26 + -26.31894506957162*t20 + 8.772981689857206*t19 + 4.386490844928603*t18)/(t3*t3*t3) + t28*t8*t8*t8*(-35.09192675942882*t8*t9 + -35.09192675942882)/(t3*t3*t3))/kT) + R*(epsilon_dispersion*kappa_association*(t1*t2*t7*rtemp10[i]*rtemp14[i]*rtemp14[i]*x[i]*(-1.5000000000000002*t10*t17*t25/((t3*t3)*(t3*t3)) + -1.0*t40) + lambda_dispersion*(rtemp10[i]*rtemp14[i]*x[i]*(0.125*t1*t2*t23*t39*t7*rtemp14[i]/(t3*t3*t3) + 0.11936620731892152*t1*t15*t17*t2*t25*t7/((t3*t3)*(t3*t3))) + lambda_dispersion*lambda_dispersion*t1*t2*t40*rtemp10[i]*rtemp14[i]*rtemp14[i]*x[i])) + R*(epsilon_dispersion*(t7*rtemp14[i]*x[i]*(-8.952465548919114e-2*t17*t22/(t16*t16*t3*t3*t3) + -0.35809862195676456*t10*t14*t36/(t16*t16*(t3*t3)*(t3*t3)))/(t1*t2*rtemp10[i]) + lambda_dispersion*(x[i]*(t37*t38*t7*rtemp14[i]/(t1*t16*t16*t2*t3*t3*t3) + 2.8496582899407503e-2*t14*t15*t36*t7/(t1*t16*t16*t2*(t3*t3)*(t3*t3)))/rtemp10[i] + lambda_dispersion*lambda_dispersion*t17*t35*rtemp14[i]*x[i]/(t1*t16*t16*t2*t3*t3*t3*rtemp10[i])))/kappa_association + R*(epsilon_dispersion*(t7*rtemp14[i]*rtemp14[i]*x[i]*(t17*((-1.5000000000000002*t20*t31 + -220.48907861546533*t5 + 0.37500000000000006*t33)/t3 + 36.74817976924422)/(t3*t3*t3) + -1.0*t34) + lambda_dispersion*(rtemp14[i]*x[i]*(6.283185307179586*t23*t42*t7*rtemp14[i]/(t3*t3*t3) + t17*t7*(0.11936620731892152*t30*t5/t3 + -2.984155182973038e-2*t30)/((t3*t3)*(t3*t3))) + lambda_dispersion*lambda_dispersion*t34*rtemp14[i]*rtemp14[i]*x[i]) + epsilon_dispersion*t7*rtemp14[i]*rtemp14[i]*x[i]*((t8*t8)*(t8*t8)*t9*(t10*(110.24453930773267*t36 + -440.97815723093066*t27)/t3 + -293.98543815395374*t28)/(t3*t3*t3) + -2.984155182973038e-2*t17*t22*(t8*t8)*(t8*t8)*t9/(t3*t3*t3))/kT) + R*(epsilon_dispersion*kappa_association*(-12.566370614359172*t1*t10*t13*t2*t25*t7*rtemp10[i]*rtemp14[i]*rtemp14[i]*rtemp14[i]*x[i]/((t3*t3)*(t3*t3)) + lambda_dispersion*t1*t13*t15*t2*t25*t7*rtemp10[i]*rtemp14[i]*rtemp14[i]*x[i]/((t3*t3)*(t3*t3))) + R*(epsilon_dispersion*(t7*rtemp14[i]*rtemp14[i]*x[i]*(-1.5*t10*t14*t17/(t16*t16*(t3*t3)*(t3*t3)) + -1.0*t24)/(t1*t2*rtemp10[i]) + lambda_dispersion*(rtemp14[i]*x[i]*(0.125*t22*t23*t7*rtemp14[i]/(t1*t16*t16*t2*t3*t3*t3) + 0.11936620731892152*t14*t15*t17*t7/(t1*t16*t16*t2*(t3*t3)*(t3*t3)))/rtemp10[i] + lambda_dispersion*lambda_dispersion*t24*rtemp14[i]*rtemp14[i]*x[i]/(t1*t2*rtemp10[i])))/kappa_association + R*(epsilon_dispersion*(t7*rtemp14[i]*rtemp14[i]*rtemp14[i]*x[i]*((-1847.164985533675*t13*t5 + 3.141592653589793*t13*t33 + -12.566370614359172*t12*t33)/t3 + 307.86083092227915*t13)/(t3*t3*t3) + lambda_dispersion*t7*rtemp14[i]*rtemp14[i]*x[i]*(-0.25*t13*t30 + t12*t30)/((t3*t3)*(t3*t3)) + epsilon_dispersion*t7*rtemp14[i]*rtemp14[i]*rtemp14[i]*x[i]*((t8*t8)*(t8*t8)*t9*(615.7216618445583*t13 + -3694.32997106735*t12)/(t3*t3*t3) + t17*(t8*t8)*(t8*t8)*t9*(1042.9786195192705*1 + -1847.164985533675*t20 + 461.79124638341875*t19 + 461.79124638341875*t18)/((t3*t3)*(t3*t3)))/kT) + R*R*(epsilon_dispersion*(-12.56637061435917*t10*t13*t14*t7*rtemp14[i]*rtemp14[i]*rtemp14[i]*x[i]/(t1*t16*t16*t2*(t3*t3)*(t3*t3)*rtemp10[i]) + lambda_dispersion*t13*t14*t15*t7*rtemp14[i]*rtemp14[i]*x[i]/(t1*t16*t16*t2*(t3*t3)*(t3*t3)*rtemp10[i]))/kappa_association + R*epsilon_dispersion*epsilon_dispersion*t10*t7*(t8*t8)*(t8*t8)*t9*(rtemp14[i]*rtemp14[i])*(rtemp14[i]*rtemp14[i])*x[i]*(3868.693299013926*t13 + -15474.773196055705*t12)/(kT*(t3*t3)*(t3*t3))))))))))))))) + epsilon_dispersion*kappa_association*lambda_dispersion*t1*t2*rtemp10[i]*x[i]*(6.332573977646111e-3*t4*t6 + 0.3183098861837907*t4*t49)/(R*R*t3*t3*t3*rtemp14[i]);
	}

	VectorXcd ktemp99(gd.NxNyNzOver2);
	ktemp99 = fft(gd, rtemp98);

	rtemp98.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp99[0] = ktemp99[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp99[i] = ktemp99[i]*exp(-2.0*lambda_dispersion*lambda_dispersion*length_scaling*length_scaling*t2*t2) + ktemp96[i]*complex(0,1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	ktemp96.resize(0); // KSpace
	VectorXd rtemp101(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT;
		const double t5 = rtemp42[i]/(t4*rtemp12[i]*rtemp12[i]);
		const double t6 = 1/rtemp10[i];
		const double t7 = 100.53096491487338*1 + 50.26548245743669*t3;
		rtemp101[i] = R*(kT*x[i]*(-1.0*t6*t7 + t5*t7) + R*R*kT*t6*x[i]*(631.6546816697189*t6 + -631.6546816697189*t5)/(kappa_association*t1*t2*t4)) + kT*kappa_association*t1*t2*x[i]*(-8.0*t3*rtemp42[i]/(t6*rtemp12[i]*rtemp12[i]) + 8.0*t3*t4)/R;
	}

	VectorXcd ktemp102(gd.NxNyNzOver2);
	ktemp102 = fft(gd, rtemp101);

	rtemp101.resize(0); // Realspace
	VectorXd rtemp103(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1/(-1.0*1 + exp(epsilon_association/kT));
		const double t2 = sqrt(0.15915494309189535*kappa_association*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R*t1) + 1) + -1.0;
		const double t3 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		const double t4 = 1/(sqrt(0.15915494309189535*kappa_association*rtemp10[i]/(R*R*t1*t3) + 1));
		rtemp103[i] = 50.26548245743669*kT*kappa_association*t4*rtemp10[i]*x[i]/(t1*t2*rtemp12[i]) + R*R*kT*t3*x[i]*(3968.8034150783765*R*R*t1*t2*t3/(kappa_association*rtemp10[i]) + -631.6546816697189*1 + -315.82734083485946*t4)/rtemp12[i];
	}

	VectorXcd ktemp104(gd.NxNyNzOver2);
	ktemp104 = fft(gd, rtemp103);

	rtemp103.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp104[0] = R*ktemp104[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp104[i] = ktemp104[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp106(gd.NxNyNz);
	rtemp106 = ifft(gd, ktemp104);

	ktemp104.resize(0); // KSpace
	VectorXd rtemp107(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1/(-1.0*1 + exp(epsilon_association/kT));
		const double t2 = sqrt(0.15915494309189535*kappa_association*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R*t1) + 1) + -1.0;
		const double t3 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		const double t4 = 1/(sqrt(0.15915494309189535*kappa_association*rtemp10[i]/(R*R*t1*t3) + 1));
		rtemp107[i] = 2.6666666666666665*kT*kappa_association*t4*rtemp10[i]*x[i]/(t1*t2*rtemp12[i]) + R*R*kT*t3*x[i]*(210.55156055657295*R*R*t1*t2*t3/(kappa_association*rtemp10[i]) + -33.510321638291124*1 + -16.755160819145562*t4)/rtemp12[i];
	}

	VectorXcd ktemp108(gd.NxNyNzOver2);
	ktemp108 = fft(gd, rtemp107);

	rtemp107.resize(0); // Realspace
	VectorXcd ktemp109(gd.NxNyNzOver2);
	ktemp109[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp109[i] = ktemp108[i]*complex(0,1)*k_i[1]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp110(gd.NxNyNz);
	rtemp110 = ifft(gd, ktemp109);

	ktemp109.resize(0); // KSpace
	VectorXcd ktemp111(gd.NxNyNzOver2);
	ktemp111[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp111[i] = ktemp108[i]*complex(0,1)*k_i[2]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp112(gd.NxNyNz);
	rtemp112 = ifft(gd, ktemp111);

	ktemp111.resize(0); // KSpace
	ktemp108[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp108[i] = ktemp108[i]*complex(0,1)*k_i[0]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp114(gd.NxNyNz);
	rtemp114 = ifft(gd, ktemp108);

	ktemp108.resize(0); // KSpace
	VectorXd rtemp115(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1/(-1.0*1 + exp(epsilon_association/kT));
		const double t2 = sqrt(0.15915494309189535*kappa_association*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R*t1) + 1) + -1.0;
		const double t3 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		const double t4 = 1/(sqrt(0.15915494309189535*kappa_association*rtemp10[i]/(R*R*t1*t3) + 1));
		rtemp115[i] = 100.53096491487338*kT*kappa_association*t4*rtemp10[i]*x[i]/(t1*t2*rtemp12[i]) + R*R*kT*t3*x[i]*(7937.606830156753*R*R*t1*t2*t3/(kappa_association*rtemp10[i]) + -1263.3093633394378*1 + -631.6546816697189*t4)/rtemp12[i];
	}

	VectorXcd ktemp116(gd.NxNyNzOver2);
	ktemp116 = fft(gd, rtemp115);

	rtemp115.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp116[0] = R*ktemp116[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp116[i] = ktemp116[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp118(gd.NxNyNz);
	rtemp118 = ifft(gd, ktemp116);

	ktemp116.resize(0); // KSpace
	VectorXd rtemp119(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp119[i] = R*kT*t4*x[i]*(631.6546816697189*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -100.53096491487338*1 + -50.26548245743669*t3)/rtemp12[i] + 8.0*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp120(gd.NxNyNzOver2);
	ktemp120 = fft(gd, rtemp119);

	rtemp119.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp120[0] = R*ktemp120[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp120[i] = ktemp120[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp122(gd.NxNyNz);
	rtemp122 = ifft(gd, ktemp120);

	ktemp120.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = log(1 + -1.0*rtemp6[i]);
		const double t2 = 1 + -1.0*rtemp6[i];
		const double t3 = t1*t2*t2 + rtemp6[i];
		const double t4 = rtemp18[i]*rtemp18[i];
		const double t5 = rtemp22[i]*rtemp22[i];
		const double t6 = rtemp2[i]*rtemp2[i];
		const double t7 = rtemp4[i]*rtemp4[i];
		const double t8 = (t3*(-1.768388256576615e-2*1/rtemp6[i] + 1.768388256576615e-2*1/t2) + 8.841941282883075e-3)/t2 + -8.841941282883075e-3*1 + 1.768388256576615e-2*rtemp25[i];
		rtemp122[i] = rtemp7[i]*(t8*rtemp118[i]*rtemp2[i]/(rtemp6[i]*rtemp6[i]) + rtemp122[i]/t2) + rtemp77[i]/t2 + (t8*rtemp106[i]*(-3.0*t7 + t6 + -3.0*t5 + -3.0*t4) + t3*(rtemp114[i]*rtemp4[i] + rtemp112[i]*rtemp22[i] + rtemp110[i]*rtemp18[i])/t2)/(t2*rtemp6[i]*rtemp6[i]) + R*(t3*rtemp2[i]*rtemp68[i]/(t2*t2*rtemp6[i]*rtemp6[i]) + kT*t3*(-0.33333333333333337*t7 + 0.33333333333333337*t6 + -0.33333333333333337*t5 + -0.33333333333333337*t4)/(t2*t2*rtemp6[i]*rtemp6[i])) + -1.0*kT*t1/R;
	}

	rtemp118.resize(0); // Realspace
	rtemp114.resize(0); // Realspace
	rtemp112.resize(0); // Realspace
	rtemp110.resize(0); // Realspace
	rtemp106.resize(0); // Realspace
	rtemp77.resize(0); // Realspace
	VectorXcd ktemp124(gd.NxNyNzOver2);
	ktemp124 = fft(gd, rtemp122);

	rtemp122.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp124[0] = R*(ktemp124[i] + 2.0*ktemp102[i]);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		const complex t2 = 1/(sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		ktemp124[i] = t2*ktemp124[i]*sin(t1) + t2*ktemp102[i]*sin(2.0*t1);
	}

	ktemp102.resize(0); // KSpace
	VectorXd rtemp126(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp126[i] = R*kT*t4*x[i]*(7937.606830156753*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -1263.3093633394378*1 + -631.6546816697189*t3)/rtemp12[i] + 100.53096491487338*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp127(gd.NxNyNzOver2);
	ktemp127 = fft(gd, rtemp126);

	rtemp126.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp127[0] = R*ktemp127[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp127[i] = ktemp127[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp129(gd.NxNyNz);
	rtemp129 = ifft(gd, ktemp127);

	ktemp127.resize(0); // KSpace
	VectorXcd ktemp130(gd.NxNyNzOver2);
	ktemp130[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp130[i] = ktemp71[i]*complex(0,1)*k_i[1];
	}

	VectorXd rtemp131(gd.NxNyNz);
	rtemp131 = ifft(gd, ktemp130);

	ktemp130.resize(0); // KSpace
	VectorXd rtemp132(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp132[i] = rtemp51[i]/(R*R);
	}

	VectorXcd ktemp133(gd.NxNyNzOver2);
	ktemp133 = fft(gd, rtemp132);

	rtemp132.resize(0); // Realspace
	VectorXcd ktemp134(gd.NxNyNzOver2);
	ktemp134[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp134[i] = ktemp133[i]*complex(0,1)*k_i[1]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp135(gd.NxNyNz);
	rtemp135 = ifft(gd, ktemp134);

	ktemp134.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp135[i] = -25.132741228718345*rtemp135[i] + rtemp131[i];
	}

	rtemp131.resize(0); // Realspace
	VectorXd rtemp137(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1);
		rtemp137[i] = 7.957747154594767e-2*kT*kappa_association*t1*rtemp10[i]*x[i]/((R*R)*(R*R)*t2*rtemp12[i]*(-1.0*1 + t2));
	}

	VectorXcd ktemp138(gd.NxNyNzOver2);
	ktemp138 = fft(gd, rtemp137);

	rtemp137.resize(0); // Realspace
	VectorXcd ktemp139(gd.NxNyNzOver2);
	ktemp139[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp139[i] = ktemp138[i]*complex(0,1)*k_i[1]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp140(gd.NxNyNz);
	rtemp140 = ifft(gd, ktemp139);

	ktemp139.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp140[i] = 50.26548245743669*rtemp140[i] + rtemp135[i];
	}

	rtemp135.resize(0); // Realspace
	VectorXd rtemp142(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t2 = 4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0;
		const double t3 = R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion);
		const double t4 = 4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + t3;
		const double t5 = -2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t6 = R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776);
		const double t7 = lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion);
		const double t8 = 2.25855*1 + t7 + t6;
		rtemp142[i] = 4.0*kT*x[i]/(R*R*rtemp12[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*(t2*(-4.1887902047863905*R*R*R*rtemp14[i]*(3.0*t5*t8/t1 + -0.5*t8) + -1.0*t5)/(t1*t1*t1) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*t5 + -4.1887902047863905*R*R*R*t2*(3.0*t4*t5/t1 + -0.5*t4))/(t1*t1*t1))/kT));
	}

	VectorXcd ktemp143(gd.NxNyNzOver2);
	ktemp143 = fft(gd, rtemp142);

	rtemp142.resize(0); // Realspace
	VectorXcd ktemp144(gd.NxNyNzOver2);
	ktemp144[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp144[i] = ktemp143[i]*complex(0,1)*k_i[1]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp145(gd.NxNyNz);
	rtemp145 = ifft(gd, ktemp144);

	ktemp144.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp145[i] = -12.566370614359172*rtemp145[i] + rtemp140[i];
	}

	rtemp140.resize(0); // Realspace
	VectorXd rtemp147(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1/(-1.0*1 + exp(epsilon_association/kT));
		const double t2 = rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT;
		rtemp147[i] = 50.26548245743669*kT*t1*x[i]*(sqrt(0.15915494309189535*kappa_association*t2*rtemp10[i]/(R*R*t1) + 1) + -1.0)/(kappa_association*t2*t2*rtemp10[i]*rtemp12[i]);
	}

	VectorXcd ktemp148(gd.NxNyNzOver2);
	ktemp148 = fft(gd, rtemp147);

	VectorXcd ktemp149(gd.NxNyNzOver2);
	ktemp149[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp149[i] = ktemp148[i]*complex(0,1)*k_i[1]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp150(gd.NxNyNz);
	rtemp150 = ifft(gd, ktemp149);

	ktemp149.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp150[i] = 6.283185307179586*rtemp150[i] + rtemp145[i];
	}

	rtemp145.resize(0); // Realspace
	VectorXcd ktemp152(gd.NxNyNzOver2);
	ktemp152[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp152[i] = ktemp133[i]*complex(0,1)*k_i[2]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp153(gd.NxNyNz);
	rtemp153 = ifft(gd, ktemp152);

	ktemp152.resize(0); // KSpace
	VectorXcd ktemp154(gd.NxNyNzOver2);
	ktemp154[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp154[i] = ktemp138[i]*complex(0,1)*k_i[2]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp155(gd.NxNyNz);
	rtemp155 = ifft(gd, ktemp154);

	ktemp154.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp155[i] = 50.26548245743669*rtemp155[i] + -25.132741228718345*rtemp153[i];
	}

	rtemp153.resize(0); // Realspace
	VectorXcd ktemp157(gd.NxNyNzOver2);
	ktemp157[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp157[i] = ktemp143[i]*complex(0,1)*k_i[2]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp158(gd.NxNyNz);
	rtemp158 = ifft(gd, ktemp157);

	ktemp157.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp158[i] = -12.566370614359172*rtemp158[i] + rtemp155[i];
	}

	rtemp155.resize(0); // Realspace
	VectorXcd ktemp160(gd.NxNyNzOver2);
	ktemp160[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp160[i] = ktemp148[i]*complex(0,1)*k_i[2]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp161(gd.NxNyNz);
	rtemp161 = ifft(gd, ktemp160);

	ktemp160.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp161[i] = 6.283185307179586*rtemp161[i] + rtemp158[i];
	}

	rtemp158.resize(0); // Realspace
	VectorXcd ktemp163(gd.NxNyNzOver2);
	ktemp163[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp163[i] = ktemp71[i]*complex(0,1)*k_i[2];
	}

	VectorXd rtemp164(gd.NxNyNz);
	rtemp164 = ifft(gd, ktemp163);

	ktemp163.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp164[i] = rtemp164[i] + rtemp161[i];
	}

	rtemp161.resize(0); // Realspace
	VectorXd rtemp166(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp166[i] = R*kT*t4*x[i]*(-421.1031211131459*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + 67.02064327658225*1 + 33.510321638291124*t3)/rtemp12[i] + -5.333333333333333*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp167(gd.NxNyNzOver2);
	ktemp167 = fft(gd, rtemp166);

	rtemp166.resize(0); // Realspace
	VectorXcd ktemp168(gd.NxNyNzOver2);
	ktemp168[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp168[i] = ktemp167[i]*complex(0,1)*k_i[1]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp169(gd.NxNyNz);
	rtemp169 = ifft(gd, ktemp168);

	ktemp168.resize(0); // KSpace
	VectorXd rtemp170(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp170[i] = R*kT*t4*x[i]*(210.55156055657295*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -33.510321638291124*1 + -16.755160819145562*t3)/rtemp12[i] + 2.6666666666666665*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp171(gd.NxNyNzOver2);
	ktemp171 = fft(gd, rtemp170);

	rtemp170.resize(0); // Realspace
	VectorXcd ktemp172(gd.NxNyNzOver2);
	ktemp172[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp172[i] = ktemp171[i]*complex(0,1)*k_i[1]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp173(gd.NxNyNz);
	rtemp173 = ifft(gd, ktemp172);

	ktemp172.resize(0); // KSpace
	VectorXd rtemp174(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp174[i] = R*kT*t4*x[i]*(421.1031211131459*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -67.02064327658225*1 + -33.510321638291124*t3)/rtemp12[i] + 5.333333333333333*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp175(gd.NxNyNzOver2);
	ktemp175 = fft(gd, rtemp174);

	rtemp174.resize(0); // Realspace
	VectorXcd ktemp176(gd.NxNyNzOver2);
	ktemp176[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp176[i] = ktemp175[i]*complex(0,1)*k_i[1]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp177(gd.NxNyNz);
	rtemp177 = ifft(gd, ktemp176);

	ktemp176.resize(0); // KSpace
	VectorXcd ktemp178(gd.NxNyNzOver2);
	ktemp178[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp178[i] = ktemp167[i]*complex(0,1)*k_i[2]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp179(gd.NxNyNz);
	rtemp179 = ifft(gd, ktemp178);

	ktemp178.resize(0); // KSpace
	VectorXcd ktemp180(gd.NxNyNzOver2);
	ktemp180[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp180[i] = ktemp171[i]*complex(0,1)*k_i[2]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp181(gd.NxNyNz);
	rtemp181 = ifft(gd, ktemp180);

	ktemp180.resize(0); // KSpace
	VectorXcd ktemp182(gd.NxNyNzOver2);
	ktemp182[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp182[i] = ktemp175[i]*complex(0,1)*k_i[2]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp183(gd.NxNyNz);
	rtemp183 = ifft(gd, ktemp182);

	ktemp182.resize(0); // KSpace
	ktemp167[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp167[i] = ktemp167[i]*complex(0,1)*k_i[0]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp185(gd.NxNyNz);
	rtemp185 = ifft(gd, ktemp167);

	ktemp167.resize(0); // KSpace
	ktemp171[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp171[i] = ktemp171[i]*complex(0,1)*k_i[0]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp187(gd.NxNyNz);
	rtemp187 = ifft(gd, ktemp171);

	ktemp171.resize(0); // KSpace
	ktemp175[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp175[i] = ktemp175[i]*complex(0,1)*k_i[0]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp189(gd.NxNyNz);
	rtemp189 = ifft(gd, ktemp175);

	ktemp175.resize(0); // KSpace
	VectorXd rtemp190(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp190[i] = R*kT*t4*x[i]*(70.18385351885766*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -11.17010721276371*1 + -5.585053606381855*t3)/rtemp12[i] + 0.8888888888888891*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp191(gd.NxNyNzOver2);
	ktemp191 = fft(gd, rtemp190);

	rtemp190.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp191[0] = R*ktemp191[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp191[i] = ktemp191[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp193(gd.NxNyNz);
	rtemp193 = ifft(gd, ktemp191);

	ktemp191.resize(0); // KSpace
	VectorXd rtemp194(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp194[i] = R*kT*t4*x[i]*(-7937.606830156753*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + 1263.3093633394378*1 + 631.6546816697189*t3)/rtemp12[i] + -100.53096491487338*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp195(gd.NxNyNzOver2);
	ktemp195 = fft(gd, rtemp194);

	rtemp194.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp195[0] = R*ktemp195[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp195[i] = ktemp195[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp197(gd.NxNyNz);
	rtemp197 = ifft(gd, ktemp195);

	ktemp195.resize(0); // KSpace
	VectorXd rtemp198(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp198[i] = R*kT*t4*x[i]*(3968.8034150783765*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -631.6546816697189*1 + -315.82734083485946*t3)/rtemp12[i] + 50.26548245743669*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp199(gd.NxNyNzOver2);
	ktemp199 = fft(gd, rtemp198);

	rtemp198.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp199[0] = R*ktemp199[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp199[i] = ktemp199[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp201(gd.NxNyNz);
	rtemp201 = ifft(gd, ktemp199);

	ktemp199.resize(0); // KSpace
	VectorXd rtemp202(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp202[i] = R*kT*t4*x[i]*(-3968.8034150783765*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + 631.6546816697189*1 + 315.82734083485946*t3)/rtemp12[i] + -50.26548245743669*kT*kappa_association*t1*t2*t3*rtemp10[i]*x[i]/(R*rtemp12[i]);
	}

	VectorXcd ktemp203(gd.NxNyNzOver2);
	ktemp203 = fft(gd, rtemp202);

	rtemp202.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp203[0] = R*ktemp203[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp203[i] = ktemp203[i]*sin(R*t1)*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1*t1)/t1;
	}

	VectorXd rtemp205(gd.NxNyNz);
	rtemp205 = ifft(gd, ktemp203);

	ktemp203.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp51[i] = rtemp51[i]/R;
	}

	VectorXcd ktemp207(gd.NxNyNzOver2);
	ktemp207 = fft(gd, rtemp51);

	rtemp51.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp207[0] = -2.0*R*ktemp207[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp207[i] = ktemp207[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp209(gd.NxNyNz);
	rtemp209 = ifft(gd, ktemp207);

	ktemp207.resize(0); // KSpace
	VectorXd rtemp210(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1);
		rtemp210[i] = 7.957747154594767e-2*kT*kappa_association*t1*rtemp10[i]*x[i]/(R*R*R*t2*rtemp12[i]*(-1.0*1 + t2));
	}

	VectorXcd ktemp211(gd.NxNyNzOver2);
	ktemp211 = fft(gd, rtemp210);

	rtemp210.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp211[0] = -2.0*R*ktemp211[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp211[i] = ktemp211[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp213(gd.NxNyNz);
	rtemp213 = ifft(gd, ktemp211);

	ktemp211.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp213[i] = -50.26548245743669*rtemp213[i] + 25.132741228718345*rtemp209[i];
	}

	rtemp209.resize(0); // Realspace
	VectorXd rtemp215(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t2 = 4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0;
		const double t3 = R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion);
		const double t4 = 4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + t3;
		const double t5 = -2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t6 = R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776);
		const double t7 = lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion);
		const double t8 = 2.25855*1 + t7 + t6;
		rtemp215[i] = 4.0*kT*x[i]/(R*rtemp12[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*(t2*(-4.1887902047863905*R*R*R*rtemp14[i]*(3.0*t5*t8/t1 + -0.5*t8) + -1.0*t5)/(t1*t1*t1) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*t5 + -4.1887902047863905*R*R*R*t2*(3.0*t4*t5/t1 + -0.5*t4))/(t1*t1*t1))/kT));
	}

	VectorXcd ktemp216(gd.NxNyNzOver2);
	ktemp216 = fft(gd, rtemp215);

	rtemp215.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp216[0] = -2.0*R*ktemp216[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp216[i] = ktemp216[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp218(gd.NxNyNz);
	rtemp218 = ifft(gd, ktemp216);

	ktemp216.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp218[i] = 12.566370614359172*rtemp218[i] + rtemp213[i];
	}

	rtemp213.resize(0); // Realspace
	VectorXd rtemp220(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1/(-1.0*1 + exp(epsilon_association/kT));
		const double t2 = rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT;
		rtemp220[i] = 50.26548245743669*R*kT*t1*x[i]*(sqrt(0.15915494309189535*kappa_association*t2*rtemp10[i]/(R*R*t1) + 1) + -1.0)/(kappa_association*t2*t2*rtemp10[i]*rtemp12[i]);
	}

	VectorXcd ktemp221(gd.NxNyNzOver2);
	ktemp221 = fft(gd, rtemp220);

	rtemp220.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp221[0] = -2.0*R*ktemp221[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp221[i] = ktemp221[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp223(gd.NxNyNz);
	rtemp223 = ifft(gd, ktemp221);

	ktemp221.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp223[i] = -6.283185307179586*rtemp223[i] + rtemp218[i];
	}

	rtemp218.resize(0); // Realspace
	ktemp133[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp133[i] = ktemp133[i]*complex(0,1)*k_i[0]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp226(gd.NxNyNz);
	rtemp226 = ifft(gd, ktemp133);

	ktemp133.resize(0); // KSpace
	ktemp138[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp138[i] = ktemp138[i]*complex(0,1)*k_i[0]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp228(gd.NxNyNz);
	rtemp228 = ifft(gd, ktemp138);

	ktemp138.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp228[i] = 50.26548245743669*rtemp228[i] + -25.132741228718345*rtemp226[i];
	}

	rtemp226.resize(0); // Realspace
	ktemp143[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp143[i] = ktemp143[i]*complex(0,1)*k_i[0]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp231(gd.NxNyNz);
	rtemp231 = ifft(gd, ktemp143);

	ktemp143.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp231[i] = -12.566370614359172*rtemp231[i] + rtemp228[i];
	}

	rtemp228.resize(0); // Realspace
	ktemp148[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp148[i] = ktemp148[i]*complex(0,1)*k_i[0]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + R*cos(t2))/t1;
	}

	VectorXd rtemp234(gd.NxNyNz);
	rtemp234 = ifft(gd, ktemp148);

	ktemp148.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp234[i] = 6.283185307179586*rtemp234[i] + rtemp231[i];
	}

	rtemp231.resize(0); // Realspace
	ktemp71[0] = 0;
	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp71[i] = ktemp71[i]*complex(0,1)*k_i[0];
	}

	VectorXd rtemp237(gd.NxNyNz);
	rtemp237 = ifft(gd, ktemp71);

	ktemp71.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		rtemp237[i] = rtemp237[i] + rtemp234[i];
	}

	rtemp234.resize(0); // Realspace
	VectorXd rtemp239(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp239[i] = kT*x[i]*(t4*(7937.606830156753*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -1263.3093633394378*1 + -631.6546816697189*t3) + 100.53096491487338*kappa_association*t1*t2*t3*rtemp10[i]/(R*R))/rtemp12[i];
	}

	VectorXcd ktemp240(gd.NxNyNzOver2);
	ktemp240 = fft(gd, rtemp239);

	rtemp239.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp240[0] = -2.0*R*ktemp240[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp240[i] = ktemp240[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp242(gd.NxNyNz);
	rtemp242 = ifft(gd, ktemp240);

	ktemp240.resize(0); // KSpace
	VectorXd rtemp243(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp243[i] = kT*x[i]*(t4*(3968.8034150783765*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + -631.6546816697189*1 + -315.82734083485946*t3) + 50.26548245743669*kappa_association*t1*t2*t3*rtemp10[i]/(R*R))/rtemp12[i];
	}

	VectorXcd ktemp244(gd.NxNyNzOver2);
	ktemp244 = fft(gd, rtemp243);

	rtemp243.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp244[0] = -2.0*R*ktemp244[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp244[i] = ktemp244[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp246(gd.NxNyNz);
	rtemp246 = ifft(gd, ktemp244);

	ktemp244.resize(0); // KSpace
	VectorXd rtemp247(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp247[i] = kT*x[i]*(t4*(-3968.8034150783765*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + 631.6546816697189*1 + 315.82734083485946*t3) + -50.26548245743669*kappa_association*t1*t2*t3*rtemp10[i]/(R*R))/rtemp12[i];
	}

	VectorXcd ktemp248(gd.NxNyNzOver2);
	ktemp248 = fft(gd, rtemp247);

	rtemp247.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp248[0] = -2.0*R*ktemp248[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp248[i] = ktemp248[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp250(gd.NxNyNz);
	rtemp250 = ifft(gd, ktemp248);

	ktemp248.resize(0); // KSpace
	VectorXd rtemp251(gd.NxNyNz);
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1) + -1.0);
		const double t3 = 1/(sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1));
		const double t4 = 1/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT);
		rtemp251[i] = kT*x[i]*(t4*(-7937.606830156753*R*R*t4/(kappa_association*t1*t2*rtemp10[i]) + 1263.3093633394378*1 + 631.6546816697189*t3) + -100.53096491487338*kappa_association*t1*t2*t3*rtemp10[i]/(R*R))/rtemp12[i];
	}

	VectorXcd ktemp252(gd.NxNyNzOver2);
	ktemp252 = fft(gd, rtemp251);

	rtemp251.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp252[0] = -2.0*R*ktemp252[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0];
		const complex t2 = R*sqrt(t1);
		ktemp252[i] = ktemp252[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*t1)*(-1.0*sin(t2)/(sqrt(t1)) + -1.0*R*cos(t2));
	}

	VectorXd rtemp254(gd.NxNyNz);
	rtemp254 = ifft(gd, ktemp252);

	ktemp252.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = -1.0*1 + exp(epsilon_association/kT);
		const double t2 = sqrt(0.15915494309189535*kappa_association*t1*rtemp10[i]*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*((4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0)/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -4.1887902047863905*R*R*R*(4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0)*(3.0*(-2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1) + -2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/((-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)*(-4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1)))/kT)/(R*R) + 1);
		const double t3 = 1/t2;
		const double t4 = -4.1887902047863905*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t5 = 4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0;
		const double t6 = R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*R*R*R*rtemp14[i]*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion);
		const double t7 = 4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + t6;
		const double t8 = -2.0943951023931953*R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855) + 1;
		const double t9 = R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776);
		const double t10 = lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion);
		const double t11 = 2.25855*1 + t9 + t10;
		rtemp147[i] = kT*x[i]*((100.53096491487338*1 + 50.26548245743669*t3)/(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*(t5*(-4.1887902047863905*R*R*R*rtemp14[i]*(3.0*t11*t8/t4 + -0.5*t11) + -1.0*t8)/(t4*t4*t4) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*t8 + -4.1887902047863905*R*R*R*t5*(3.0*t7*t8/t4 + -0.5*t7))/(t4*t4*t4))/kT) + -8.0*kappa_association*t1*t3*rtemp10[i]/(R*R*(-1.0*1 + t2)))/(R*R*rtemp12[i]) + -12.566370614359172*rtemp147[i];
	}

	VectorXcd ktemp256(gd.NxNyNzOver2);
	ktemp256 = fft(gd, rtemp147);

	rtemp147.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp256[0] = R*(-2.0*ktemp79[i] + ktemp256[i]);
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = sin(R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t2 = 1/(sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]));
		const complex t3 = t1*t2;
		ktemp256[i] = ktemp79[i]*(-1.0*R*cos(R/t2) + -1.0*t3) + t3*ktemp256[i];
	}

	ktemp79.resize(0); // KSpace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp256[0] = kT*ktemp1[i]/(R*R) + ktemp256[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp256[i] = ktemp256[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0])) + kT*ktemp1[i]/(R*R);
	}

	ktemp1.resize(0); // KSpace
	VectorXd rtemp259(gd.NxNyNz);
	rtemp259 = ifft(gd, ktemp256);

	ktemp256.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = 1 + -1.0*rtemp6[i];
		const double t2 = -7.957747154594767e-2*rtemp4[i]*rtemp4[i] + 7.957747154594767e-2*rtemp2[i]*rtemp2[i] + -7.957747154594767e-2*rtemp22[i]*rtemp22[i] + -7.957747154594767e-2*rtemp18[i]*rtemp18[i];
		const double t3 = rtemp18[i]*rtemp18[i];
		const double t4 = rtemp22[i]*rtemp22[i];
		const double t5 = rtemp2[i]*rtemp2[i];
		const double t6 = rtemp4[i]*rtemp4[i];
		const double t7 = -3.0*t6 + t5 + -3.0*t4 + -3.0*t3;
		const double t8 = t1*t1*log(t1) + rtemp6[i];
		const double t9 = -0.22222222222222224*1/rtemp6[i] + 0.22222222222222224*1/t1;
		const double t10 = log(t1);
		const double t11 = t8/(t1*rtemp6[i]) + -1.0*rtemp25[i];
		const double t12 = (7.957747154594767e-2*t8*t9 + 8.841941282883075e-3)/t1;
		rtemp259[i] = ((0.3333333333333333*t2*rtemp246[i] + 0.3333333333333333*t2*(t8*rtemp254[i]/t1 + rtemp250[i])/t1)/t1 + 0.3333333333333333*t11*t2*rtemp242[i]/t1)/(rtemp6[i]*rtemp6[i]) + rtemp2[i]*((t7*(rtemp197[i]*(-8.841941282883075e-3*1 + 1.768388256576615e-2*rtemp25[i] + t12) + t8*rtemp193[i]/(t1*rtemp6[i]))/(t1*rtemp6[i]) + 7.957747154594767e-2*t7*t9*(t10*rtemp197[i] + rtemp205[i])/t1 + t7*rtemp201[i]*(-8.841941282883075e-3*1 + 1.768388256576615e-2*rtemp25[i] + 7.957747154594767e-2*t9 + 2.0*t12)/(t1*t1) + t7*rtemp193[i]*(t8/(t1*t1) + 1)/(t1*t1))/(rtemp6[i]*rtemp6[i]) + rtemp4[i]*(((t8*rtemp189[i]/t1 + rtemp187[i])/t1 + rtemp50[i])/t1 + t11*rtemp185[i]/t1)/(rtemp6[i]*rtemp6[i]) + rtemp22[i]*(((t8*rtemp183[i]/t1 + rtemp181[i])/t1 + rtemp92[i])/t1 + t11*rtemp179[i]/t1)/(rtemp6[i]*rtemp6[i]) + rtemp18[i]*(((t8*rtemp177[i]/t1 + rtemp173[i])/t1 + rtemp85[i])/t1 + t11*rtemp169[i]/t1)/(rtemp6[i]*rtemp6[i]) + rtemp223[i]/(t1*t1)) + rtemp259[i]/t1 + rtemp237[i]*rtemp4[i]/(t1*t1) + rtemp164[i]*rtemp22[i]/(t1*t1) + rtemp150[i]*rtemp18[i]/(t1*t1) + kT*rtemp7[i]*(-0.22222222222222224*t10*t7 + (t7*t8*t9 + 0.11111111111111112*t7)/t1 + -0.11111111111111112*t7)/(rtemp6[i]*rtemp6[i]) + (1.5*rtemp2[i]*rtemp68[i]/(t1*t1) + t2*(rtemp129[i]/t1 + 12.566370614359172*kT)/(t1*t1))/R;
	}

	rtemp254.resize(0); // Realspace
	rtemp250.resize(0); // Realspace
	rtemp246.resize(0); // Realspace
	rtemp242.resize(0); // Realspace
	rtemp237.resize(0); // Realspace
	rtemp223.resize(0); // Realspace
	rtemp205.resize(0); // Realspace
	rtemp201.resize(0); // Realspace
	rtemp197.resize(0); // Realspace
	rtemp193.resize(0); // Realspace
	rtemp189.resize(0); // Realspace
	rtemp187.resize(0); // Realspace
	rtemp185.resize(0); // Realspace
	rtemp183.resize(0); // Realspace
	rtemp181.resize(0); // Realspace
	rtemp179.resize(0); // Realspace
	rtemp177.resize(0); // Realspace
	rtemp173.resize(0); // Realspace
	rtemp169.resize(0); // Realspace
	rtemp164.resize(0); // Realspace
	rtemp150.resize(0); // Realspace
	rtemp129.resize(0); // Realspace
	rtemp92.resize(0); // Realspace
	rtemp85.resize(0); // Realspace
	rtemp68.resize(0); // Realspace
	rtemp50.resize(0); // Realspace
	rtemp25.resize(0); // Realspace
	rtemp22.resize(0); // Realspace
	rtemp18.resize(0); // Realspace
	rtemp7.resize(0); // Realspace
	rtemp6.resize(0); // Realspace
	rtemp4.resize(0); // Realspace
	rtemp2.resize(0); // Realspace
	VectorXcd ktemp261(gd.NxNyNzOver2);
	ktemp261 = fft(gd, rtemp259);

	rtemp259.resize(0); // Realspace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp261[0] = 0.33333333333333337*R*R*R*ktemp261[i] + ktemp124[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		const complex t1 = R*sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0]);
		ktemp261[i] = ktemp261[i]*(-1.0*t1*cos(t1) + sin(t1))/(sqrt(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0])*(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0])) + ktemp124[i];
	}

	ktemp124.resize(0); // KSpace
		{
		const int i = 0;
		const Reciprocal k_i = Reciprocal(0,0,0);
		ktemp261[0] = ktemp99[i] + ktemp261[i];
	}

	for (int i=1; i<gd.NxNyNzOver2; i++) {
		const int z = i % gd.NzOver2;
		const int n = (i-z)/gd.NzOver2;
		const int y = n % gd.Ny;
		const int xa = (n-y)/gd.Ny;
		const RelativeReciprocal rvec((xa>gd.Nx/2) ? xa - gd.Nx : xa, (y>gd.Ny/2) ? y - gd.Ny : y, z);
		const Reciprocal k_i = gd.Lat.toReciprocal(rvec);
		ktemp261[i] = ktemp261[i]*exp(-6.0*pow(gd.dvolume, 0.6666666666666666)*(k_i[2]*k_i[2] + k_i[1]*k_i[1] + k_i[0]*k_i[0])) + ktemp99[i];
	}

	ktemp99.resize(0); // KSpace
	VectorXd rtemp264(gd.NxNyNz);
	rtemp264 = ifft(gd, ktemp261);

	ktemp261.resize(0); // KSpace
	for (int i=0; i<gd.NxNyNz; i++) {
		const double t1 = R*R*R*rtemp14[i]*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(4.1887902047863905*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 17.54596337971441*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -2.803431620357388) + 2.25855);
		const double t2 = 1 + -4.1887902047863905*t1;
		const double t3 = -4.1887902047863905*R*R*R*rtemp14[i]*(3.0*(1 + -2.0943951023931953*t1)*(lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + 2.25855)/t2 + -0.5*lambda_dispersion*(-1.50349*1 + 0.249434*lambda_dispersion) + -0.5*R*R*R*rtemp14[i]*(8.377580409572781*lambda_dispersion*(1.40049*1 + -0.827739*lambda_dispersion) + 52.63789013914324*R*R*R*rtemp14[i]*(lambda_dispersion*(-15.0427*1 + 5.30827*lambda_dispersion) + 10.1576) + -5.606863240714776) + -1.129275) + -1.0*1 + 2.0943951023931953*t1;
		const double t4 = R*R*R*rtemp14[i];
		const double t5 = 4.1887902047863905*t4;
		const double t6 = 4.0*lambda_dispersion*lambda_dispersion*lambda_dispersion + -4.0;
		const double t7 = 12.566370614359172*R*R*(sqrt(0.15915494309189535*kappa_association*rtemp10[i]*(-1.0*1 + exp(epsilon_association/kT))*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*(t3*t6/(t2*t2*t2) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(1 + -2.0943951023931953*t1) + -4.1887902047863905*R*R*R*t6*(-2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + 3.0*(1 + -2.0943951023931953*t1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*t4*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/t2 + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*t4*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/(t2*t2*t2))/kT)/(R*R) + 1) + -1.0)/(kappa_association*rtemp10[i]*(-1.0*1 + exp(epsilon_association/kT))*(rtemp42[i]/rtemp12[i] + 0.25*epsilon_dispersion*(t3*t6/(t2*t2*t2) + -7.957747154594767e-2*lambda_dispersion*(-50.26548245743669*lambda_dispersion*lambda_dispersion*(1 + -2.0943951023931953*t1) + -4.1887902047863905*R*R*R*t6*(-2.0943951023931953*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + 3.0*(1 + -2.0943951023931953*t1)*(4.1887902047863905*rtemp14[i]*(-1.50349*1 + 0.498868*lambda_dispersion) + R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*t4*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion))/t2 + -0.5*R*R*R*rtemp14[i]*rtemp14[i]*(73.49635953848843*t4*(-15.0427*1 + 10.61654*lambda_dispersion) + 24.572946253656237*1 + -29.046956363922856*lambda_dispersion)))/(t2*t2*t2))/kT));
		(*outgrad)[i] = ingrad[i]*(kT*(2.0*1 + log(2.6464769766182683e-6*x[i]/(sqrt(kT)*kT)) + 4.0*log(t7) + -2.0*t7) + epsilon_dispersion*t4*t6*(2.0943951023931953*epsilon_dispersion*t3*((1 + -1.0*t5)*(1 + -1.0*t5))*((1 + -1.0*t5)*(1 + -1.0*t5))/(kT*(4.0*t5*t5 + 1 + 4.0*t5)) + -4.1887902047863905*1 + 8.772981689857206*t1)/(t2*t2*t2) + rtemp264[i] + mu) + (*outgrad)[i];
	}

	rtemp264.resize(0); // Realspace
	rtemp42.resize(0); // Realspace
	rtemp14.resize(0); // Realspace
	rtemp12.resize(0); // Realspace
	rtemp10.resize(0); // Realspace
	// 119 Fourier transform used.
	// 44

}

void print_summary(const char *prefix, double energy, std::string name) const {
	if (name != "") printf("%s%25s =", prefix, name.c_str());
	else printf("%s%25s =", prefix, "UNKNOWN");
	print_double("", energy);
	printf("\n%s%25s =", prefix, "FSAFT");
	print_double("", FSAFT);
	printf("\n%s%25s =", prefix, "Fassoc");
	print_double("", Fassoc);
	printf("\n%s%25s =", prefix, "Fdisp");
	print_double("", Fdisp);
	printf("\n%s%25s =", prefix, "Fideal");
	print_double("", Fideal);
	printf("\n%s%25s =", prefix, "a1integrated");
	print_double("", a1integrated);
	printf("\n%s%25s =", prefix, "a2integrated");
	print_double("", a2integrated);
	printf("\n%s%25s =", prefix, "dV");
	print_double("", dV);
	printf("\n%s%25s =", prefix, "dr");
	print_double("", dr);
	printf("\n%s%25s =", prefix, "volume");
	print_double("", volume);
	printf("\n%s%25s =", prefix, "whitebear");
	print_double("", whitebear);
	printf("\n");
}

private:
	double R;
	double epsilon_association;
	double kappa_association;
	double epsilon_dispersion;
	double lambda_dispersion;
	double length_scaling;
	double mu;
	mutable double FSAFT;
	mutable double Fassoc;
	mutable double Fdisp;
	mutable double Fideal;
	mutable double a1integrated;
	mutable double a2integrated;
	mutable double dV;
	mutable double dr;
	mutable double volume;
	mutable double whitebear;
	// TODO: add declaration of spherical fourier transform data here
	mutable double oldkT;}; // End of WaterSaft_type class
	// Total 139 Fourier transform used.
	// peak memory used: 44



Functional WaterSaft(double R, double epsilon_association, double kappa_association, double epsilon_dispersion, double lambda_dispersion, double length_scaling, double mu) {
	return Functional(new WaterSaft_type(R, epsilon_association, kappa_association, epsilon_dispersion, lambda_dispersion, length_scaling, mu), "WaterSaft");
}

