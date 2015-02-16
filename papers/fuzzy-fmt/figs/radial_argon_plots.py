#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
from pylab import *

sigma = 0.3506 #nm
sigma_over_R=2**(5/6)

figure()
data = loadtxt('figs/YarnellArgon85K.dat')
n =0.02125 # Angstrom^-3
nsig_3 = n*(sigma*10)**3
plot(data[:,0]/(sigma*10),data[:,1])
data_mc1=loadtxt('figs/mc_testp_wca-0.8389-0.7100.dat.gradial')
plot((1/sigma_over_R)*data_mc1[:,0],data_mc1[:,1])
data_dft1 = loadtxt('figs/radial-lj-0.7100-0.84.dat')
plot(data_dft1[:,0]/sigma_over_R,data_dft1[:,1])
xlabel(r'$r/\sigma$')
savefig('figs/Argon-vapor_pressure-85K.pdf')

figure()
data2 = loadtxt('figs/EggertArgon0.6GPaX.dat')
n = 24.23 #atoms/nm^3
nsig_3 = n*sigma**3
plot(data2[:,0]/sigma,data2[:,1])
data_mc2 = loadtxt('figs/mc_testp_wca-0.9570-2.4800.dat.gradial')
plot((1/sigma_over_R)*data_mc2[:,0], data_mc2[:,1])
data_dft2 = loadtxt('figs/radial-lj-2.4800-0.96.dat')
plot(data_dft2[:,0]/sigma_over_R,data_dft2[:,1])
xlabel(r'$r/\sigma$')
savefig('figs/Argon-0_6GPa-RT.pdf')

figure()
data3 = loadtxt('figs/EggertArgon1.1GPaRAW.dat')
n = 27.74 #atoms/nm^3

nsig_3 = n*sigma**3
plot(data3[:,0]/sigma,(data3[:,1]))
data_mc3 = loadtxt('figs/mc_testp_wca-1.0950-2.4800.dat.gradial')
plot((1/sigma_over_R)*data_mc3[:,0], data_mc3[:,1])
data_dft3 = loadtxt('figs/radial-lj-2.4800-1.09.dat')
plot(data_dft3[:,0]/sigma_over_R,data_dft3[:,1])
xlabel(r'$r/\sigma$')
savefig('figs/Argon-1_1GPa-RT.pdf')

figure()
data4 = loadtxt('figs/Mikolaj-X.dat')
#n = 27.74 #atoms/nm^3
#nsig_3 = n*sigma**3
data_mc4 = loadtxt('figs/mc_testp_wca-0.5844-1.2350.dat.gradial')
plot(data4[:,0]/(sigma*10),data4[:,1])
plot((1/sigma_over_R)*data_mc4[:,0], data_mc4[:,1])
data_dft4 = loadtxt('figs/radial-lj-1.2350-0.58.dat')
plot(data_dft4[:,0]/sigma_over_R,data_dft4[:,1]*1.7)
xlabel(r'$r/\sigma$')


savefig('figs/Argon-9_919MPa-148K.pdf')
