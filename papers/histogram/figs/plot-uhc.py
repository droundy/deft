#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

import styles

if len(sys.argv) not in [5,6]:
    print 'useage: %s ww ff N versions show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.3]

# note: speficic HC should be independent of N, but we have to choose one
N = float(sys.argv[3])
#arg N = range(5,21)

versions = eval(sys.argv[4])
#arg versions = [["nw","wang_landau","simple_flat","optimized_ensemble","kT0.4","kT0.5","tmmc","oetmmc"]]

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat" % (ww, ff, N, version, data) for version in versions for data in ["E","lnw"]]

reference_method = "wang_landau"

max_T = 2
T_bins = 1e3
dT = max_T/T_bins
T_range = numpy.arange(dT,max_T,dT)
min_T = 0 # we will adjust this

fig_u = plt.figure('u')
plt.title('Specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

fig_hc = plt.figure('hc')
plt.title('Specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

# make dictionaries which we can index by method name
U = {} # internal energy
CV = {} # heat capacity

for version in versions:

    with open("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version)) as file:
        for line in file:
            if("min_T" in line):
                this_min_T = float(line.split()[-1])
                if this_min_T > min_T:
                    min_T = this_min_T
                break

    # energy histogram file; indexed by [-energy,counts]
    e_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, version), ndmin=2)
    # weight histogram file; indexed by [-energy,ln(weight)]
    lnw_hist = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-%s-lnw.dat" % (ww, ff, N, version), ndmin=2)

    energy = -e_hist[:,0] # array of energies
    lnw = lnw_hist[e_hist[:,0].astype(int),1] # look up the lnw for each actual energy
    ln_dos = numpy.log(e_hist[:,1]) - lnw

    Z = numpy.zeros(len(T_range)) # partition function
    U[version] = numpy.zeros(len(T_range)) # internal energy
    CV[version] = numpy.zeros(len(T_range)) # heat capacity
    for i in range(len(T_range)):
        ln_dos_boltz = ln_dos - energy/T_range[i]
        dos_boltz = numpy.exp(ln_dos_boltz - ln_dos_boltz.max())
        Z[i] = sum(dos_boltz)
        U[version][i] = sum(energy*dos_boltz)/Z[i]
        CV[version][i] = sum((energy/T_range[i])**2*dos_boltz)/Z[i] - \
                         (sum(energy/T_range[i]*dos_boltz)/Z[i])**2

    plt.figure('u')
    plt.plot(T_range,U[version]/N,styles.plot(version),label=styles.title(version))

    plt.figure('hc')
    plt.plot(T_range,CV[version]/N,styles.plot(version),label=styles.title(version))

plt.figure('u')
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u.pdf" % (ww*100, ff*100, N))

plt.figure('hc')
plt.ylim(0)
plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc.pdf" % (ww*100, ff*100, N))

min_T_i = int(min_T/max_T*T_bins)
error_data = open("figs/error-table-ww%02.0f-ff%02.0f-%i.dat" % (ww*100, ff*100, N), "w")
error_data.write("# method u_error cv_error\n")
error_data.write("# min_T: %g\n" % min_T)
for version in versions:
    u_error = max(abs(U[version][min_T_i:] - U[reference_method][min_T_i:]))/N
    cv_error = max(abs(CV[version][min_T_i:] - CV[reference_method][min_T_i:]))/N
    error_data.write("%s %g %g\n" % (version, u_error, cv_error))
error_data.close()

if 'show' in sys.argv:
    plt.show()

