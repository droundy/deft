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
eref, lndosref = readnew.e_lndos(ref)
maxref = readnew.max_entropy_state(ref)

r = glob('data/%s-movie/*lndos.dat' % (filebase))

iterations = numpy.zeros(len(r))
maxentropystate = numpy.zeros(len(r))
minimportantenergy = numpy.zeros(len(r))
erroratenergy = numpy.zeros(len(r))

goodenough = 0.1
goodenoughenergy = numpy.zeros(len(r))

for i,f in enumerate(sorted(r)):
    e,lndos = readnew.e_lndos(f)
    iterations[i] = readnew.iterations(f)
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
plt.plot(iterations, erroratenergy, 'r.')
plt.xlabel('# iterations')
plt.ylabel('error')

plt.figure()
plt.plot(iterations, goodenoughenergy, 'r.')
plt.xlabel('# iterations')
plt.ylabel('energy we are converged to')
plt.title('Convergence at %g%% level' % (goodenough*100))
plt.show()

















