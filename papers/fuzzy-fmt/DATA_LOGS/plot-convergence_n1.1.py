#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

#gw=0.10 Phi_1
plt.figure()
print "n1.1 gw 0.10 Phi_1"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.10_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.10_fv_0_vs_Phi_1.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,2]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_1 vs dx at n1.1 gw0.10"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_1$ a.k.a. Phi 1')
plt.legend()
#plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.10  Phi_2
plt.figure()
print "n1.1 gw 0.10 Phi_2"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.10_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.10_fv_0_vs_Phi_2.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,3]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_2 vs dx at n1.1 gw0.10"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_2$ a.k.a. Phi 2')
plt.legend()
#plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.10 Phi_3
plt.figure()
print "n1.1 gw 0.10 Phi_3"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.10_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.10_fv_0_vs_Phi_3.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,4]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_3 vs dx at n1.1 gw0.10"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_3$ a.k.a. Phi 3')
plt.legend()
#plt.legend(loc='best')

plt.savefig(plot_name)

#---------------------------

#gw=0.13 Phi_1
plt.figure()
print "n1.1 gw 0.13 Phi_1"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.13_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.13_fv_0_vs_Phi_1.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,2]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_1 vs dx at n1.1 gw0.13"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_1$ a.k.a. Phi 1')
plt.legend()
#plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.13 Phi_2
plt.figure()
print "n1.1 gw 0.13 Phi_2"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.13_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.13_fv_0_vs_Phi_2.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,3]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_2 vs dx at n1.1 gw0.13"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_2$ a.k.a. Phi 2')
plt.legend()
#plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.13 Phi_3
plt.figure()
print "n1.1 gw 0.13 Phi_3"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.13_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.13_fv_0_vs_Phi_3.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,4]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_3 vs dx at n1.1 gw0.13"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_3$ a.k.a. Phi 3')
plt.legend()
#plt.legend(loc='best')

plt.savefig(plot_name)

#---------------------------
#gw=0.15  Phi_1
plt.figure()
print "n1.1 gw 0.15 Phi_1"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.15_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.15_fv_0_vs_Phi_1.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,2]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_1 vs dx at n1.1 gw0.15"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_1$ a.k.a. Phi 1')
plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.15  Phi_2
plt.figure()
print "n1.1 gw 0.15 Phi_2"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.15_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.15_fv_0_vs_Phi_2.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,3]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_2 vs dx at n1.1 gw0.15"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_2$ a.k.a. Phi 2')
plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.15  Phi_3
plt.figure()
print "n1.1 gw 0.15 Phi_3"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.15_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.15_fv_0_vs_Phi_3.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,4]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_3 vs dx at n1.1 gw0.15"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_3$ a.k.a. Phi 3')
plt.legend(loc='best')

plt.savefig(plot_name)

#---------------------------
#gw=0.18  Phi_1
plt.figure()
print "n1.1 gw 0.18 Phi_1"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.18_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.18_fv_0_vs_Phi_1.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,2]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_1 vs dx at n1.1 gw0.18"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_1$ a.k.a. Phi 1')
plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.18  Phi_2
plt.figure()
print "n1.1 gw 0.18 Phi_2"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.18_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.18_fv_0_vs_Phi_2.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,3]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_2 vs dx at n1.1 gw0.18"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_2$ a.k.a. Phi 2')
plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.18  Phi_3
plt.figure()
print "n1.1 gw 0.18 Phi_3"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.18_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.18_fv_0_vs_Phi_3.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,4]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_3 vs dx at n1.1 gw0.18"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_3$ a.k.a. Phi 3')
plt.legend(loc='best')

plt.savefig(plot_name)

#---------------------------
#gw=0.20  Phi_1
plt.figure()
print "n1.1 gw 0.20 Phi_1"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.20_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.20_fv_0_vs_Phi_1.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,2]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_1 vs dx at n1.1 gw0.20"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_1$ a.k.a. Phi 1')
plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.20  Phi_2
plt.figure()
print "n1.1 gw 0.20 Phi_2"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.20_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.20_fv_0_vs_Phi_2.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,3]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_2 vs dx at n1.1 gw0.20"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_2$ a.k.a. Phi 2')
plt.legend(loc='best')

plt.savefig(plot_name)

#gw=0.20  Phi_3
plt.figure()
print "n1.1 gw 0.20 Phi_3"
data = np.loadtxt('DATA_LOGS/n1.1_T2_gw_0.20_fv_0.dat')
plot_name="DATA_LOGS/PLOTS/Plots_n1.1/n1.1_T2_gw_0.20_fv_0_vs_Phi_3.png"

dx = data[:,0]
mcerror = data[:,1]
phi1 = data[:,4]

mce_values = set()
for e in mcerror:
    mce_values.add(e)
#print(mce_values)

for e in sorted(mce_values):
    plt.plot(dx[mcerror==e], phi1[mcerror==e], 'x', label='mcerror %s' % e)

plot_title="Phi_3 vs dx at n1.1 gw0.20"
plt.title(plot_title)
plt.xlabel(r'$\Delta x$')
plt.ylabel(r'$\Phi_1$ a.k.a. Phi 3')
plt.legend(loc='best')

plt.savefig(plot_name)

#use plt.show to display all figures with zoom option
#plt.show()
