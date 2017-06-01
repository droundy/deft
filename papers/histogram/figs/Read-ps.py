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

methods = [ '-tmi3', '-tmi2', '-tmi', '-toe', '-tmmc']

#----------------------------------------------------------------------#
nom = len(methods)

T = numpy.zeros(nom,dtype = object)
u = numpy.zeros(nom,dtype = object)
cv = numpy.zeros(nom,dtype = object)
s = numpy.zeros(nom,dtype = object)
minT = numpy.zeros(nom,dtype = object)

energy = numpy.zeros(nom,dtype = object)
lndos = numpy.zeros(nom,dtype = object)
ps = numpy.zeros(nom,dtype = object)

fref = 'data/lv/ww%.2f-ff%.2f-%gx%g-tmi3-dos.dat' % (ww,ff,lenx,lenyz)
energy_ref,lndos_ref = readnew.e_lndos(fref)
max_entropy_state = readnew.max_entropy_state(fref)

for i in range(len(methods)):
    try:
        method = methods[i]
        fbase = 'data/lv/ww%.2f-ff%.2f-%gx%g%s-movie/%s' % (ww,ff,lenx,lenyz,method,number)
        T[i],u[i],cv[i],s[i],minT[i] = readnew.T_u_cv_s_minT(fbase)
        energy[i],lndos[i],ps[i] = readnew.e_lndos_ps(fbase)
        energy[i] = energy[i][:len(lndos_ref)]
        lndos[i] = lndos[i][:len(lndos_ref)]
        ps[i] = ps[i][:len(lndos_ref)]

    except:
        pass

#----------------------------------------------------------------------#

Error_DOS = numpy.absolute(lndos[0] - lndos[1])
#plt.semilogx(ps1[0], Error_DOS, label = '$\mid tmi_3 - tmi_2 \mid$')
error0 = abs(1-numpy.exp(lndos[0]-lndos[0][max_entropy_state]
                         -lndos_ref+lndos_ref[max_entropy_state])[ps[0] > 0])

error1 = abs(1-numpy.exp(lndos[1]-lndos[1][max_entropy_state]
                         -lndos_ref+lndos_ref[max_entropy_state])[ps[1] > 0])

plt.loglog(ps[0][ps[0] > 0], error0,
             'ro',linewidth = 2, label = '$tmi_3$')
plt.loglog(ps[1][ps[1] > 0], error1,
             'bx',linewidth = 2, label = '$tmi_2$')
plt.loglog(ps[1][ps[1] > 0], 3/np.sqrt(ps[1][ps[1] > 0]), 'k:',
            label=r'$\frac{3}{\sqrt{ps}}$')
#plt.ylim((-50,-20))

plt.xlabel(r'Pessimistic Samples')
plt.ylabel('$(D-D_{ref})/D_{ref}$')
plt.title('$\mathcal{D}(\epsilon)$ for $\lambda=%g$, $\eta=%g$' % (ww, ff))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g%s-Error.pdf'
            % (ww,ff,lenx,lenyz,'-tmi3'))

plt.figure()
#plt.semilogx(ps1[0], Error_DOS, label = '$\mid tmi_3 - tmi_2 \mid$')
plt.plot(energy[0], lndos[0],'r',linewidth = 2, label = '$tmi_3$')
plt.plot(energy[1], lndos[1],'b',linewidth = 2, label = '$tmi_2$')

plt.xlabel(r'energy')
plt.ylabel('$\mathcal{D}(\epsilon)$')
plt.title('$\mathcal{D}(\epsilon)$ for $\lambda=%g$, $\eta=%g$' % (ww, ff))
plt.legend(loc='best')
plt.tight_layout(pad=0.2)

plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g%s-vs-energy.pdf'
            % (ww,ff,lenx,lenyz,'-tmi3'))
plt.show()
