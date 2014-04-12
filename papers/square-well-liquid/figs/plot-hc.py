#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy

if len(sys.argv) != 2:
    print 'useage: %s ww' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

# todo: make inputs?
ktemin = 0
ktemax = 30
dkte = float(ktemax-ktemin) / 1000
kte = numpy.arange(ktemin+dkte,ktemax,dkte)

# interntal energy relative to well depth
def Ue(kte,counts,DS):
    output = numpy.zeros(len(kte))
    for i in range(len(output)):
        output[i] = sum(counts*DS*numpy.exp(-counts/kte[i])) \
          / sum(DS*numpy.exp(-counts/kte[i]))
    return output

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
fmt.set_powerlimits((-2,3))
fmt.set_scientific(True)
ax.yaxis.set_major_formatter(fmt)
ax.xaxis.set_major_formatter(fmt)

plt.title('Heat capacity for well width %g' % (ww))

# todo: figure out how to specify which ff and N to use
# input: "data/periodic-ww%04.2f-ff0.10-N200-nw-dos.dat" % (ww)
# input: "data/periodic-ww%04.2f-ff0.20-N200-nw-dos.dat" % (ww)
# input: "data/periodic-ww%04.2f-ff0.30-N200-nw-dos.dat" % (ww)
for ff in [ 0.1, 0.2, 0.3 ]:
    for N in [ 200 ]:
        data = numpy.loadtxt(
            "data/periodic-ww%04.2f-ff%04.2f-N%i-nw-dos.dat" % (ww, ff, N),
            ndmin=2)
        counts = data[:,0][::-1]
        counts -= min(counts)
        DS = data[:,1]
        cv = (Ue(kte+dkte/2,counts,DS)
              - Ue(kte-dkte/2,counts,DS)) / dkte / N
        plt.plot(kte,cv,label=r'$\eta=%g$ and $N=%i$' % (ff, N))

plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-hc.pdf" % (ww*100))
