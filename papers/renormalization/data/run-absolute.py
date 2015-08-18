#!/usr/bin/python

from __future__ import division
from math import pi       # REALLY don't need all of math
import os, numpy, sys

os.system("fac ../../../free-energy-monte-carlo")

R=1

i= 0
# RG recursion level

ww=1.3
#arg ww = [1.3, 1.5, 2.0, 3.0]

L=5
#arg L = [5.0]

Ns = list(xrange(0,4))
#arg Ns = [range(2,10)]

Overwrite = False
if '-O' in sys.argv:  # Check for overwrite flag in arguments
    Overwrite = True

def free_energy_over_kT_to_ff(free_energy):
    # bisection?
    n = 0
    while n < 1000: #better value?
        ff = np.sqrt(free_energy) # This is not TRUE!! Or approximate!
        # C-S: (4*ff - 3*ff**2)/(1-ff)**2 = F  
        if 
    return ff

for N in Ns:
    dirname = 'scrunched-ww%4.2f-L%04.2f/i%01d/N%03d/absolute' % (ww,L,i,N)
    os.system('mkdir -p '+dirname)
    if Overwrite:
        os.system('rm -r '+dirname)

    success_ratios = []
    all_total_checks = []
    all_valid_checks = []
    steps = 20 # Need a better value for this
    step_size = 0.05 # This too
    steps = 20
    sim_iters = 1000000
    sc_period = int(max(10, 1*N*N/10))

    ffs = np.zeros(steps+1)
    ffs[0] = (4/3.0)*pi*R**3/L**3
    for j in xrange(steps):
        if j!= 0 and j % 4 == 0:
            step_size = step_size * 0.5
        ffs[j+1] = ffs[j] + step_size

    # approximate_free_energies = np.arange(np.log(2), 100, np.log(2))
    # ffs = np.zeros_like(approximate_free_energies)
    # for i in xrange(len(ffs)):
    #     ffs[i] = free_energy_over_kT_to_ff(approximate_free_energies[i])

    for j in xrange(len(ffs)-1):
        filename = '%05d' % (j)

        if not os.path.isfile(dirname+'/'+filename+'.dat'):
            #cmd = 'srun -J %s' % filename
            cmd = ' ../../../free-energy-monte-carlo'
            cmd += ' --ff %g' % ffs[j]
            cmd += ' --sc_period %d' % sc_period
            cmd += ' --iterations %d' % sim_iters
            cmd += ' --filename %s' % filename
            cmd += ' --data_dir %s' % dirname
            cmd += ' --ff_small %g' % ffs[j+1]
            cmd += ' --N %d' % N
            cmd += ' > %s/%s.out 2>&1 &' % (dirname, filename)
            print("Running with command: %s" % cmd)
            os.system(cmd)
        else:
            print("You're trying to overwrite files, use the flag -O in order to do so.")
