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

# note: RD should be independent of ff and N, but we have to choose values for now
ff = float(sys.argv[2])
#arg ff = [0.3]

N = float(sys.argv[3])
#arg N = [200]

kTs = eval(sys.argv[4])
#arg kTs = [[0.1, 1, 2]]

# FIXME: make this an input?
max_radius = 7

# radial distribution function at fixed temperature
def g(gs,kte,counts,dos):
    output = numpy.zeros(len(gs[0,:]))
    Z = sum(dos*numpy.exp(-counts/kte))
    for i in range(len(output)):
        output[i] += sum(gs[:,i]*dos*numpy.exp(-counts/kte)) / Z
    return output

plt.title('Radial distribution function for $\lambda=%g$, $\eta=%g$, and $N=%i$'
          % (ww, ff, N))

# FIXME: would it be better to use data from simulations with fixed kT?

g_data = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-nw-g.dat" % (ww, ff, N), ndmin=2)
g_counts = g_data[:,0][::-1]
gs = g_data[:,1:]

e_hist = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-nw-E.dat" % (ww, ff, N), ndmin=2)
lnw_hist = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-nw-lnw.dat" % (ww, ff, N), ndmin=2)
counts = e_hist[:,0][::-1]
use_counts = [ i for i in range(len(counts))
               if counts[i] in g_counts ]
counts = counts[use_counts]
counts -= min(counts)
lnw_mean = lnw_hist[:,1].mean()
dos_data = e_hist[:,1]*numpy.exp(-(lnw_hist[:,1] - lnw_mean))
dos = dos_data[use_counts]
dos /= sum(dos)

with open("data/periodic-ww%04.2f-ff%04.2f-N%i-nw-g.dat"
          % (ww, ff, N),'r') as stream:
    first_line = stream.readline().split(' ')
for i in range(len(first_line)):
    if 'de_g' in first_line[i]:
        de_g = float(first_line[i+1])
        break
radius = (numpy.array(range(0,len(gs[0,:])))+0.5) * de_g/2

for kT in kTs:
    plt.plot(radius,g(gs,kT,counts,dos),'.',
             label='$kT/\epsilon=%g$' %(kT))

plt.axvline(1,color='k',linestyle=':')
plt.axvline(ww,color='k',linestyle=':')
plt.xlim(0,max_radius/2)
plt.xlabel('$r/\\sigma$')
plt.ylabel('$g(r)$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-rdf.pdf" % (ww*100, ff*100, N))
