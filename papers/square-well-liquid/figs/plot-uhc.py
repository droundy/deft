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
kTmax = 10
dkT = float(kTmax-kTmin) / 1000
kT_range = numpy.arange(kTmin+dkT,kTmax,dkT)

# make figure with axes labeled using scientific notation
def sci_fig():
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-2,3))
    fmt.set_scientific(True)
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_formatter(fmt)
    return fig,ax

#############################
### internal energy plots ###
#############################

# internal energy per ball as a function of kT
def u(kT_array,data):
    energy = -data[:,0]
    DS = data[:,1]
    u_out = numpy.zeros(len(kT_array))
    for i in range(len(u_out)):
        u_out[i] = sum(energy*DS*numpy.exp(-(energy-min(energy))/kT_array[i])) \
          / sum(DS*numpy.exp(-(energy-min(energy))/kT_array[i]))
    return u_out/N

fig, ax = sci_fig()
plt.title('Specific internal energy for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

# plot curves from simulations for all temperatures
data = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-flat-dos.dat" % (ww, ff, N),
                     ndmin=2)
plt.plot(kT_range,u(kT_range,data),label='flat histogram')

# plot curves from simulations at fixed temperature
# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-kT%g-dos.dat" % (ww, ff, N, kT) for kT in kTs]
for kT in kTs:
    data = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-kT%g-dos.dat" % (ww, ff, N, kT),
        ndmin=2)
    plt.plot(kT_range,u(kT_range,data),label='$kT=%g\epsilon$ sim.' % kT)

plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-u.pdf" % (ww*100, ff*100, N))
plt.close()

###########################
### heat capacity plots ###
###########################

# specific heat capacity as a function of kT
def cv(data):
    return (u(kT_range+dkT/2,data) - u(kT_range-dkT/2,data)) / dkT

fig, ax = sci_fig()
plt.title('Specific heat capacity for $\lambda=%g$, $\eta=%g$, and $N=%i$' % (ww, ff, N))

# plot curves from simulations for all temperatures
data = numpy.loadtxt("data/periodic-ww%04.2f-ff%04.2f-N%i-flat-dos.dat" % (ww, ff, N),
                     ndmin=2)
plt.plot(kT_range,cv(data),label='flat histogram')

# plot curves from simulations at fixed temperature
# input: ["data/periodic-ww%04.2f-ff%04.2f-N%i-kT%g-dos.dat" % (ww, ff, N, kT) for kT in kTs]
for kT in kTs:
    data = numpy.loadtxt(
        "data/periodic-ww%04.2f-ff%04.2f-N%i-kT%g-dos.dat" % (ww, ff, N, kT),
        ndmin=2)
    plt.plot(kT_range,cv(data),label='$kT=%g\epsilon$ sim.' % kT)

plt.xlabel('$kT/\epsilon$')
plt.ylabel('$C_V/Nk$')
plt.legend(loc='best')
plt.tight_layout(pad=0.1)
plt.savefig("figs/periodic-ww%02.0f-ff%02.0f-N%i-hc.pdf" % (ww*100, ff*100, N))
plt.close()
