#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import readandcompute

if len(sys.argv) not in [4,5]:
    print 'useage: %s ww ff Ns show' % sys.argv[0]
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3, 1.5, 2.0, 3.0]

ff = float(sys.argv[2])
#arg ff = [0.3]

# note: speficic HC should be independent of N, but we have to choose one
Ns = eval(sys.argv[3])
#arg Ns = [[500,1372,4000]]

# make figure with axes labeled using scientific notation
def sci_fig(handle):
    fig = plt.figure(handle)
    ax = fig.add_subplot(1,1,1)
    fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-2,3))
    fmt.set_scientific(True)
    ax.yaxis.set_major_formatter(fmt)
    ax.xaxis.set_major_formatter(fmt)
    return fig, ax

fig_u, ax_u = sci_fig('u')
plt.title('Specific internal energy for $\lambda=%g$, $\eta=%g$' % (ww, ff))

for N in Ns:
    T, U, CV, S, minT = readandcompute.T_u_cv_s_minT("data/mc/ww%.2f-ff%.2f-N%d" % (ww, ff, N))
    plt.plot(T,U/N,'-', label='$N=%d' % N)

plt.xlabel('$kT/\epsilon$')
plt.ylabel('$U/N\epsilon$')
plt.legend(loc='best')
plt.tight_layout(pad=0.2)
plt.savefig("figs/energy-ww%02.0f-ff%02.0f.pdf" % (ww*100, ff*100))

plt.show()
