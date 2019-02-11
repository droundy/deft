from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys, glob, matplotlib
import colors
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')

filename = sys.argv[1]

plt.figure(figsize=(5,4))
print('starting!')
try:
    for wl in glob.glob("data/gamma/%s/wl*.txt" % filename):
        print "in wl"
        print 'vanilla_wang_landau'+ wl[len("data/gamma/%s/wl" % filename):-4]
        wlmoves, wlfactor = np.loadtxt(wl, dtype = float, unpack = True)
        data = np.loadtxt(wl)
        moves = data[:,0]
        factor = data[:,1]
        if (data[0,0] == 'wl_factor'): # using c++ data!
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
        plt.ylim(ymin=1e-10, ymax=1e1)
        plt.xlim(xmin=1e0, xmax=2e12)

except:
    pass

print(glob.glob("data/gamma/%s/sad*.dat" % filename))
try:
    for sad in glob.glob("data/gamma/%s/sad*.dat" % filename):
        data = np.loadtxt(sad)
        ts = data[:,0]
        avg_gamma = data[:,1]
        sadname = sad.split('/')[-1].split('.')[0]
    
    
        plt.ylim(ymin=1e-10, ymax=1e1)
        plt.xlim(xmin=1e0, xmax=2e12)
        #print(avg_gamma)
        #print(ts)
        colors.loglog(ts, avg_gamma,sadname)
        if data.shape[1] > 2:
            max_avg_gamma = data[:,2]
            min_avg_gamma = data[:,3]
            plt.fill_between(ts, min_avg_gamma, max_avg_gamma,
                            edgecolor='none', linewidth=0,
                            color=colors.color('sad'),
                            alpha=0.1, zorder=-51)
except:
    raise

def gamma_sa(t,t0):
    return t0/np.maximum(t, t0)

t0s = ['1e3','1e4','1e5','1e6','1e7']



for t0 in t0s:
    colors.loglog(ts,gamma_sa(ts, float(t0)),'samc-%s-%s' %(t0,filename.replace('n','')))
    plt.xlabel(r'$\textrm{Moves}$')
    plt.ylabel(r'$\gamma_{t}$')
    colors.legend()
plt.tight_layout()
plt.savefig('figs/gamma-%s.pdf' % filename.replace('.','_'))
if 'noshow' not in sys.argv:
    plt.show()

