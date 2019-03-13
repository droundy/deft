from __future__ import division, print_function
import numpy as np
import sys, os, matplotlib
from collections import OrderedDict

matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font', family='serif')
matplotlib.rcParams['figure.figsize'] = (5, 7)

if 'noshow' in sys.argv:
        matplotlib.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
import colors

if os.path.exists('../data'):
    os.chdir('..')

# Arguments for entropy portion
energy = int(sys.argv[1])
filebase = sys.argv[2]
transcale = sys.argv[3] 

tex_filebase = filebase.replace('.','_') # latex objects to extra "." characters

methods = ['-sad3','-vanilla_wang_landau']
if 'allmethods' not in sys.argv:
    methods = ['-sad3','-wl-50','-wl-inv-t-50','-wl-inv-t-256','-sad-50',
               '-sad-256','-wl-256']
    if transcale == 'slow':
        methods = ['-sad3-slow','-wl-50-slow','-wl-inv-t-50-slow','-sad-50-slow',
                   '-vanilla_wang_landau-slow']
    if transcale == 'fast':
        methods = ['-sad3-fast','-vanilla_wang_landau-fast']

# For SAMC compatibility with LVMC
lvextra1 = glob('data/comparison/%s-samc*' % filebase)
split3 = [i.split('%s-'%filebase, 1)[-1] for i in lvextra1]
split4 = [i for i in split3 if i[-3:-1] != '-s']

for meth in split4:
    if "default" in transcale:
        if meth[-3:] != '-tm' and "slow" not in meth and "fast" not in meth:
            methods.append('-%s' % meth)
    if "slow" in transcale: 
        if meth[-3:] != '-tm' and "slow" in meth:
            methods.append('-%s' % meth)
    if "fast" in transcale:
        if meth[-3:] != '-tm' and "fast" in meth:
            methods.append('-%s' % meth)

best_ever_max = 1e100
max_time = 0
min_error = 1e200
print('methods are', methods)
for method in methods:
    print('trying method', method)
    try:
        dirname = 'data/comparison/%s%s/' % (filebase,method)
        if not os.path.exists(dirname) or os.listdir(dirname) == []:
                continue

        if energy > 0:
                Nrt_at_energy, erroratenergy = np.loadtxt(dirname + 'energy-%s.txt' % energy, delimiter = '\t', unpack = True)
        data = np.loadtxt(dirname + 'errors.txt', delimiter = '\t', unpack = True)
        iterations = data[0]
        errorinentropy = data[1]
        maxerror = data[2]
        best_ever_max = min(best_ever_max, maxerror.min())

        if not os.path.exists('figs/s000'):
                os.makedirs('figs/s000')

        if filebase.startswith('s000'):
                N = filebase.split('-N')[-1]
                ff = filebase.split('-ff')[-1].split('-N')[0]
                ff = float(ff)
                # Get N directly from title.
                moves = iterations * float(N)
        max_time = max(max_time, moves.max())

        if type(iterations) is not np.float64:
                plt.figure('errorinentropy')
                ax1 = plt.subplot(2, 1, 1)
                if data.shape[0] > 4:
                    minmean = data[3]
                    maxmean = data[4]
                    plt.fill_between(moves, minmean, maxmean,
                                     edgecolor='none', linewidth=0,
                                     color=colors.color(method[1:]),
                                     alpha=0.1, zorder=-51)
                my_S_error = errorinentropy[0:len(iterations)]
                min_error = min(min_error, my_S_error[my_S_error > 0].min())

                colors.loglog(moves, my_S_error, method = method[1:])

                # make these tick labels invisible
                plt.setp(ax1.get_xticklabels(), visible=False)
                #plt.xlabel(r'$\textrm{Moves}$')
                plt.ylabel(r'$\textrm{Average Entropy Error}$')
                #plt.title('Average Entropy Error at Each Iteration, %s' %filebase)
                if "default" in transcale:
                    plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.05$' % (int(N), ff))
                elif "slow" in transcale:
                    plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.005$' % (int(N), ff))
                elif "fast" in transcale:
                    plt.title(r'$N=%d$, $\eta = %g$ $\delta_0 = 0.5$' % (int(N), ff))
                #colors.legend()

    except:
        raise
 
plt.figure('errorinentropy')
moves = np.array([1e3, 1e12])
#colors.loglog(moves, min_error*np.sqrt(moves.max())/np.sqrt(moves), method = r'1/sqrt(t)')
for i in np.arange(-8, 19, 1.0):
    colors.loglog(moves, 10**i/np.sqrt(0.1*moves), method = r'1/sqrt(t)')
plt.xlim(moves[0], moves[1])
if filebase == 's000/periodic-ww1.30-ff0.30-N50':
    plt.ylim(1e-3, 1e4)
    if "slow" in transcale:
        plt.ylim(1e-2, 1e5)
elif filebase == 's000/periodic-ww1.30-ff0.30-N500':
    plt.ylim(1e-1, 1e3)
elif filebase == 's000/periodic-ww1.50-ff0.17-N256':
    plt.ylim(1e-3, 1e3)
#colors.legend()
#plt.tight_layout()
#print('filename', 'figs/%s-entropy-error-%s.pdf' % (tex_filebase,transcale))
#plt.savefig('figs/%s-entropy-error-%s.pdf' % (tex_filebase,transcale))

#----------------------------------------------------------------------#
ax2 = plt.subplot(2, 1, 2,sharex=ax1)
ax2.xaxis.set_ticks_position('both')
if "default" in transcale:
    wl_globstring = 'wl*[!slow].txt'
    sad_globstring = 'sad*[!slow].dat'
    plt.text(10**3, 10**12.3, '(a)',fontsize=14)
    plt.text(10**3, 10**0.3, '(b)',fontsize=14)
elif "slow" in transcale:
    wl_globstring = 'wl*slow.txt'
    sad_globstring = 'sad*slow.dat'
    plt.text(10**3, 10**12.3, '(c)',fontsize=14)
    plt.text(10**3, 10**0.3, '(d)',fontsize=14)

try:
    for wl in glob("data/gamma/n%s/%s" % (N,wl_globstring)):
        print('in wl')
        print('vanilla_wang_landau'+ wl[len("data/gamma/n%s/wl" % N):-4])
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
                      'wl'
                         + wl[len("data/gamma/n%s/wl" % N):-4])
        plt.ylim(ymin=1e-10, ymax=1e1)
        plt.xlim(xmin=1e2, xmax=1e12)

except:
    pass

try:
    for sad in glob("data/gamma/n%s/%s" % (N,sad_globstring)):
        data = np.loadtxt(sad)
        #print(data.shape)
        if data.shape[1] > 2:
                print('You must not be using the rust code!?!')
                exit(1)
        else: #data.shape[0] == 2: # we are using parse-yaml-out.py
            ts = data[:,0]
            gamma = data[:,1]
            sadname = sad.split('/')[-1].split('.')[0]
    
    
        plt.ylim(ymin=1e-10, ymax=1e1)
        plt.xlim(xmin=1e3, xmax=1e12)
        print(sadname)
        colors.loglog(ts, gamma,sadname)
except:
    pass

def gamma_sa(t,t0):
    return t0/np.maximum(t, t0)

t0s = ['1e3','1e4','1e5','1e6','1e7']



for t0 in t0s:
    colors.loglog(ts,gamma_sa(ts, float(t0)),'samc-%s-%s' %(t0,N))
    plt.xlabel(r'$t$')
    plt.ylabel(r'$\gamma_{t}$')
    colors.legend()
plt.tight_layout()
plt.subplots_adjust(hspace=0.0)

print('filename', 'figs/%s-gamma-and-entropy-error-%s.pdf' % (tex_filebase,transcale))
plt.savefig('figs/%s-gamma-and-entropy-error-%s.pdf' % (tex_filebase,transcale))
if 'noshow' not in sys.argv:
    plt.show()
