#!/usr/bin/python2
import matplotlib, sys
import numpy

if len(sys.argv) != 6:
    print 'useage: %s ww ff Ns T method' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.3]

Ns = eval(sys.argv[3])
#arg Ns = [range(5,31)]

T = float(sys.argv[4])
#arg T = [0.1,0.2]

method = sys.argv[5]
#arg method = ['tmmc-golden']

# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, method) for N in Ns]

min_e = {}
for N in Ns:
    with open("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat" % (ww, ff, N, method)) as file:
        for line in file:
            if 'min_T' in line:
                golden_min_T = float(line.split()[-1])
                break

    if T < golden_min_T:
        print("We do not have golden data for N = %i and T < %g" %(N,golden_min_T))
        exit(1)

    e_hist = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-E.dat"
                           % (ww, ff, N, method))
    lnw_hist = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-%s-lnw.dat"
                             % (ww, ff, N, method))
    lnw = lnw_hist[e_hist[:,0].astype(int),1]
    ln_dos = numpy.log(e_hist[:,1]) - lnw

    for i in range(len(ln_dos)-1):
        if abs(ln_dos[i+1]-ln_dos[i]) >= 1/T:
            min_e[N] = int(e_hist[i,0])
            break
    if N not in min_e.keys():
        min_e[N] = int(e_hist[-1,0])

results = open("figs/min_important_energies-ww%02.0f-ff%02.0f-T%g.dat"
               % (ww*100, ff*100, T), "w")
results.write("# N, min_e\n")
for N in Ns:
    results.write("%i, %i\n"%(N,min_e[N]))
