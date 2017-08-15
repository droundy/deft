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
methods = [ '-tmi3', '-tmi2', '-tmi', '-toe', '-toe3','-tmmc', '-vanilla_wang_landau']
ref = "data/" + reference
maxref = readnew.max_entropy_state(ref)
goodenough = 0.1

colors = {
    '-tmi3': 'r',
    '-tmi2': 'k',
    '-tmi': 'y',
    '-toe': 'm',
    '-toe3': 'g',
    '-tmmc': 'b',
    '-toe2': 'c',
    '-vanilla_wang_landau': 'b'}
try:
    eref, lndosref, Nrt_ref = readnew.e_lndos_ps(ref)
except:
    eref, lndosref = readnew.e_lndos(ref)

for method in methods:
    try: 
        r = glob('data/%s%s-movie/*lndos.dat' % (filebase,method))
        if len(r)==0:
            continue
        iterations = numpy.zeros(len(r))
        Nrt_at_energy = numpy.zeros(len(r))
        maxentropystate = numpy.zeros(len(r))
        minimportantenergy = numpy.zeros(len(r))
        erroratenergy = numpy.zeros(len(r))
        goodenoughenergy = numpy.zeros(len(r))

        for i,f in enumerate(sorted(r)):
            e,lndos,Nrt = readnew.e_lndos_ps(f)

            iterations[i] = readnew.iterations(f)
            Nrt_at_energy[i] = Nrt[energy]
            maxentropystate[i] = readnew.max_entropy_state(f)
            minimportantenergy[i] = readnew.min_important_energy(f)

            doserror = lndos - lndos[maxref] - lndosref + lndosref[maxref]
            erroratenergy[i] = doserror[energy]
        print e
        newiterations = []
        newnrt = []
        newmax = []
        newmin = []
        i = 0
        while iterations[i] < max(iterations):
            newiterations = iterations[:i+1]
            newmin = minimportantenergy[:i+1]
            newmax = maxentropystate[:i+1]
            i+=1

        i = 0
        while Nrt_at_energy[i] < max(Nrt_at_energy):
            newnrt = Nrt_at_energy[:i+1]
            i+=1
            
        # The following finds out what energy we are converged to at the
        # "goodenough" level.
        if any(numpy.logical_and(abs(doserror) > goodenough, e < -maxref)):
            goodenoughenergy[i] = max(e[numpy.logical_and(abs(doserror) > goodenough, e < -maxref)])
        else:
            goodenoughenergy[i] = -minimportantenergy[i]
        
        plt.figure('error-at-energy-iterations')
        plt.plot(newiterations[newnrt > 0], erroratenergy[newnrt > 0], '%s-' % colors[method], label = method[1:])
        plt.title('Error at energy %g' % energy)
        plt.xlabel('# iterations')
        plt.ylabel('error')
        plt.legend(loc = 'best')

        plt.figure('round-trips-at-energy')
        plt.plot(newiterations, newnrt, '%s-' % colors[method], label = method[1:])
        plt.title('Roundy Trips at energy %g' % energy)
        plt.xlabel('# iterations')
        plt.ylabel('Roundy Trips')
        plt.legend(loc = 'best')
        
        plt.figure('error-at-energy-round-trips')
        plt.plot(newnrt[newnrt > 0], erroratenergy[newnrt > 0], '%s-' % colors[method], label = method[1:])
        plt.title('Error at energy %g' % energy)
        plt.xlabel('Roundy Trips')
        plt.ylabel('Error')
        plt.legend(loc = 'best')
        
        plt.figure('convergence')
        plt.plot(newiterations, goodenoughenergy[0:len(newiterations)], '%s-' % colors[method], label = method[1:])
        plt.xlabel('# iterations')
        plt.ylabel('energy we are converged to')
        plt.title('Convergence at %g%% level' % (goodenough*100))
        plt.legend(loc = 'best')

    except:
        print 'I had trouble with', method
        raise
plt.show()

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


















