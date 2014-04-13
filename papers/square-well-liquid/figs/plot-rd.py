#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

if len(sys.argv) != 3:
    print 'useage: %s ww N' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

# FIXME(?): RD be independent of N, but we have to choose one for now
N = float(sys.argv[2])
#arg N = [200, 1000]

# FIXME: figure out how to specify which kTs to use
kTs = [0.1, 1, 2]

# FIXME: it shouldn't matter what filling fraction we use,
#    so I'm choosing one for now.
#    We can use data from several for better statistics probably.
ff = 0.3

# FIXME: make this an inputs?
max_radius = 7

# radial distribution function at fixed temperature
def g(gs,kte,counts,DS):
    output = numpy.zeros(len(gs[0,:]))
    Z = sum(DS*numpy.exp(-counts/kte))
    for i in range(len(output)):
        output[i] += sum(gs[:,i]*DS*numpy.exp(-counts/kte)) / Z
    return output

plt.title('Heat capacity for well width %g' % (ww))

# FIXME: see fixmes about kTs and ff
# FIXME: would it be better to use data from simulations with fixed kT?
# input: "data/periodic-ww%04.2f-ff0.30-N%i-nw-g.dat" % (ww, N)
# input: "data/periodic-ww%04.2f-ff0.30-N%i-nw-dos.dat" % (ww, N)
for kT in kTs:
    g_data = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-nw-g.dat" % (ww, ff, N),
        ndmin=2)
    g_counts = g_data[:,0][::-1]
    gs = g_data[:,1:]

    dos_data = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-nw-dos.dat" % (ww, ff, N),
        ndmin=2)
    counts = dos_data[:,0][::-1]
    use_counts = [ i for i in range(len(counts))
                   if counts[i] in g_counts ]
    counts = counts[use_counts]
    counts -= min(counts)
    DS = dos_data[use_counts,1]
    DS /= sum(DS)

    with open("data/periodic-ww%04.2f-ff%04.2f-N%i-nw-g.dat"
              % (ww, ff, N),'r') as stream:
        first_line = stream.readline().split(' ')
    for i in range(len(first_line)):
        if 'de_g' in first_line[i]:
            de_g = float(first_line[i+1])
            break
    radius = (numpy.array(range(0,len(gs[0,:])))+0.5) * de_g/2

    plt.plot(radius,g(gs,kT,counts,DS),'.',
             label='$kT/\epsilon=%g$' %(kT))

plt.axvline(1,color='k',linestyle=':')
plt.axvline(ww,color='k',linestyle=':')
plt.xlim(0,max_radius/2)
plt.xlabel('$r/\\sigma$')
plt.ylabel('$g(r)$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-N%i-rd.pdf" % (ww*100, N))
