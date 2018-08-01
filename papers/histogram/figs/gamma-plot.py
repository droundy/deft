from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys, glob
import colors

from matplotlib.colors import LightSource

densitycolormap = plt.cm.jet
densityinterpolation = 'bilinear'
densityshadedflag = True
densitybarflag = True
gridflag = True

filename = sys.argv[1]
Tmin = float(sys.argv[2])

try:
    for wl in glob.glob("data/gamma/%s/wl*.txt" % filename):
        wlmoves, wlfactor = np.loadtxt(wl, dtype = float, unpack = True)
        moves = np.zeros(len(wlmoves)*2+2)
        factor = np.zeros_like(moves)
        factor[0] = 1
        moves[0] = 1
        for i in range(len(wlmoves)):
            moves[2*i+1] = wlmoves[i]
            moves[2*i+2] = wlmoves[i]
            factor[2*i+1] = wlfactor[i]*2
            factor[2*i+2] = wlfactor[i]
        colors.loglog(moves, factor,
                      'vanilla_wang_landau'
                         + wl[len("data/gamma/%s/wl" % filename):-4])
        plt.ylim(ymin=1e-10)

except:
    pass


for sad in glob.glob("data/gamma/%s/sad*.dat" % filename):
    data = np.loadtxt(sad)
    energies_found = data[:,0]
    time = data[:,1]
    ehi = data[:,2]
    elo = data[:,3]
    ts = np.exp(np.linspace(0, np.log(max(time)*5), 1000))
    gamma = np.zeros_like(ts)
    print sad, time
    for j in range(len(time)):
        for i in range(len(gamma)):
            if ts[i] > time[j]:
                t = ts[i]
                tL = time[j]
                NE = energies_found[j]
                gamma[i] = abs(elo[j]-ehi[j])/(3*Tmin*t)*(
                  NE**2 + NE*t + (t/tL-1)*t)/(NE**2 + t + (t/tL-1)*t)
    sadname = sad.split('/')[-1].split('.')[0]


    colors.loglog(ts, gamma,sadname)

def gamma_sa(t,t0):
    return t0/np.maximum(t, t0)

t0s = [1e3,1e4,1e5,1e6]
for t0 in t0s:
    colors.loglog(ts,gamma_sa(ts, t0),'samc-%g' %t0)
    plt.xlabel('Moves')
    plt.ylabel('Gamma')
    colors.legend()
plt.savefig('figs/gamma-%s.pdf' % filename.replace('.','_'))
if 'noshow' not in sys.argv:
    plt.show()

