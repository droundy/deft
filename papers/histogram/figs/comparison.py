from __future__ import division
import numpy, sys, os
import matplotlib.pyplot as plt
import readnew
from glob import glob

if os.path.exists('../data'):
    os.chdir('..')

energy = int(sys.argv[1])
reference = sys.argv[2]
filebase = sys.argv[3]

ref = "data/" + reference
try:
    eref, lndosref, Nrt_ref = readnew.e_lndos_ps(ref)
except:
    eref, lndosref = readnew.e_lndos(ref)
maxref = readnew.max_entropy_state(ref)

r = glob('data/%s-movie/*lndos.dat' % (filebase))

iterations = numpy.zeros(len(r))
Nrt_at_energy = numpy.zeros(len(r))
maxentropystate = numpy.zeros(len(r))
minimportantenergy = numpy.zeros(len(r))
erroratenergy = numpy.zeros(len(r))

goodenough = 0.1
goodenoughenergy = numpy.zeros(len(r))

for i,f in enumerate(sorted(r)):
    try:
        e,lndos,Nrt = readnew.e_lndos_ps(f)
    except:
        continue
    e,lndos = readnew.e_lndos(f)
    iterations[i] = readnew.iterations(f)
    Nrt_at_energy[i] = Nrt[energy]
    maxentropystate[i] = readnew.max_entropy_state(f)
    minimportantenergy[i] = readnew.min_important_energy(f)
    doserror = lndos - lndos[maxref] - lndosref + lndosref[maxref]
    erroratenergy[i] = doserror[energy]

    # The following finds out what energy we are converged to at the
    # "goodenough" level.
    if any(numpy.logical_and(abs(doserror) > goodenough, e < -maxref)):
        goodenoughenergy[i] = max(e[numpy.logical_and(abs(doserror) > goodenough, e < -maxref)])
    else:
        goodenoughenergy[i] = -minimportantenergy[i]
print minimportantenergy[-1]



#~ plt.loglog(ps[ps > 0], iterations[ps > 0],
             #~ 'ro',linewidth = 2,label = r'30x4 tmi3')
#~ plt.loglog(psref[ps > 0], iterations[ps > 0],
             #~ 'bx',linewidth = 2,label = r'30x4 ref')

#~ plt.xlabel(r'$iterations$')
#~ plt.ylabel(r'$N_R(\epsilon)$')
#~ #plt.title(r'$D(\epsilon)$ for $\lambda=%g$, $\eta=%g$' % (ww, ff))
#~ plt.legend(loc='best')
#~ plt.tight_layout(pad=0.2)

#plt.savefig('figs/sticky-wall-ww%.2f-ff%.2f-%gx%g%s-Error.pdf'
            #% (ww,ff,lenx,lenyz,'-tmi3'))

plt.figure()
plt.plot(iterations[Nrt_at_energy > 0], erroratenergy[Nrt_at_energy > 0], 'r.')
plt.title('Error at energy %g' % energy)
plt.xlabel('# iterations')
plt.ylabel('error')

plt.figure()
plt.plot(iterations, Nrt_at_energy, 'r.')
plt.title('Roundy Trips at energy %g' % energy)
plt.xlabel('# iterations')
plt.ylabel('Roundy Trips')

plt.figure()
plt.plot(Nrt_at_energy[Nrt_at_energy > 0], erroratenergy[Nrt_at_energy > 0], 'r.')
plt.title('Error at energy %g' % energy)
plt.xlabel('Roundy Trips')
plt.ylabel('Error')

plt.figure()
plt.plot(iterations, goodenoughenergy, 'r.')
plt.xlabel('# iterations')
plt.ylabel('energy we are converged to')
plt.title('Convergence at %g%% level' % (goodenough*100))
plt.show()

















