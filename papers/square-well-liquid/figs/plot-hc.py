#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

if len(sys.argv) != 5:
    print 'useage: %s ww ff N kTs' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.3]

# note: speficic HC should be independent of N, but we have to choose one
N = float(sys.argv[3])
#arg N = [200]

kTs = eval(sys.argv[4])
#arg kTs = [[0.1, 1, 2]]

# FIXME: make these inputs?
kTmin = 0
kTmax = 30
dkT = float(kTmax-kTmin) / 1000
kT_range = numpy.arange(kTmin+dkT,kTmax,dkT)

# interntal energy relative to well depth
def U(kT_array,counts,DS):
    output = numpy.zeros(len(kT_array))
    for i in range(len(output)):
        output[i] = sum(counts*DS*numpy.exp(-counts/kT_array[i])) \
          / sum(DS*numpy.exp(-counts/kT_array[i]))
    return output

# label axes using scientific notation
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((-2,3))
fmt.set_scientific(True)
ax.yaxis.set_major_formatter(fmt)
ax.xaxis.set_major_formatter(fmt)

# plot curves using density of state data from simulations for all temperatures
plt.title('Heat capacity for $ww=%g$, $ff=%g$, and $N=%i$' % (ww, ff, N))
data = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-flat-dos.dat" % (ww, ff, N),
                     ndmin=2)
counts = data[:,0][::-1]
counts -= min(counts)
DS = data[:,1]
cv = (U(kT_range+dkT/2,counts,DS)
      - U(kT_range-dkT/2,counts,DS)) / dkT / N
plt.plot(kT_range,cv,label=r'$\eta=%g$ and $N=%i$' % (ff, N))

# FIXME: plot points using density of state data from fixed kT simulations
#   it looks like I need to run two simulations for each point kT, namely
#     one at kT-dkT/2 and kT+dkT/2, as heat capacity is given by dU/dkT

plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc.pdf" % (ww*100, ff*100, N))
