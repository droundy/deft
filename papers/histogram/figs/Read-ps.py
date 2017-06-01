#!/usr/bin/python2
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import numpy as np

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readnew

if len(sys.argv) < 5:
    print("Usage: python {} 1.3 0.22 100 10 tmi".format(sys.argv[0]))
    exit(1)

ww = float(sys.argv[1])
#arg ww = [1.3]
ff = float(sys.argv[2])
#arg ff = [0.1, 0.2, 0.3]
lenx = float(sys.argv[3])
#arg lenx = [50, 80, 100]
lenyz = float(sys.argv[4])
#arg lenyz = [10]
number = sys.argv[5]
#arg number = [000100]

all_methods = [ '-tmi3', '-tmi2', '-tmi', '-toe', '-tmmc']
methods = []
prettymethods = []

T = []
u = []

energy = []
lndos = []
ps = []
hist = []

fref = 'data/lv/ww%.2f-ff%.2f-%gx%g-tmi3-dos.dat' % (ww,ff,lenx,lenyz)
energy_ref,lndos_ref = readnew.e_lndos(fref)
max_entropy_state = readnew.max_entropy_state(fref)

for method in all_methods:
    try:
        fbase = 'data/lv/ww%.2f-ff%.2f-%gx%g%s-movie/%s' % (ww,ff,lenx,lenyz,method,number)
        myT,myu,mycv,mys,myminT = readnew.T_u_cv_s_minT(fbase)
        myenergy,mylndos,myps = readnew.e_lndos_ps(fbase)
        _, myhist = readnew.e_hist(fbase)

    except:
        continue
    energy.append(myenergy[:len(lndos_ref)])
    lndos.append(mylndos[:len(lndos_ref)])
    ps.append(myps[:len(lndos_ref)])
    newhist = np.zeros_like(lndos_ref)
    newhist[:len(myhist)] = myhist
    hist.append(newhist)
    methods.append(method)
    prettymethods.append(method[1:])
    T.append(myT)
    u.append(myu)


error0 = abs(1-numpy.exp(lndos[0]-lndos[0][max_entropy_state]
                         -lndos_ref+lndos_ref[max_entropy_state])[ps[0] > 0])

error1 = abs(1-numpy.exp(lndos[1]-lndos[1][max_entropy_state]
                         -lndos_ref+lndos_ref[max_entropy_state])[ps[1] > 0])

plt.loglog(ps[0][ps[0] > 0], error0,
             'ro',linewidth = 2, label = prettymethods[0])
plt.loglog(ps[1][ps[1] > 0], error1,
             'bx',linewidth = 2, label = prettymethods[1])
plt.loglog(ps[1][ps[1] > 0], 3/np.sqrt(ps[1][ps[1] > 0]), 'k:',
            label=r'$\frac{3}{\sqrt{N_R}}$')
#plt.ylim((-50,-20))

plt.xlabel(r'$N_R(\epsilon)$')
plt.ylabel(r'$\frac{D(\epsilon)-D_{\textrm{ref}}(\epsilon)}{D_{\textrm{ref}}(\epsilon)}$')
plt.title(r'$D(\epsilon)$ for $\lambda=%g$, $\eta=%g$' % (ww, ff))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g%s-Error.pdf'
            % (ww,ff,lenx,lenyz,'-tmi3'))

plt.figure()
plt.plot(energy[0], lndos[0],'r',linewidth = 2, label = prettymethods[0])
plt.plot(energy[1], lndos[1],'b',linewidth = 2, label = prettymethods[1])

plt.xlabel(r'energy')
plt.ylabel('$D(\epsilon)$')
plt.title('$D(\epsilon)$ for $\lambda=%g$, $\eta=%g$' % (ww, ff))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g%s-vs-energy.pdf'
            % (ww,ff,lenx,lenyz,'-tmi3'))

plt.figure()
#print ps[0].shape, hist[0].shape
#print hist[0]
plt.loglog(ps[0][ps[0] > 0], hist[0][ps[0]>0],
             'ro',linewidth = 2, label = prettymethods[0])
plt.loglog(ps[1][ps[1] > 0], hist[1][ps[1]>0],
             'bx',linewidth = 2, label = prettymethods[1])

plt.xlabel(r'$N_R(\epsilon)$')
plt.ylabel(r'$H(\epsilon)$')
plt.title('$H(\epsilon)$ for $\lambda=%g$, $\eta=%g$' % (ww, ff))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g%s-hist-vs-NR.pdf'
            % (ww,ff,lenx,lenyz,'-tmi3'))
plt.show()
