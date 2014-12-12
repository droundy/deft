from __future__ import division
from os import system
import os
from math import pi

figsdir = 'papers/fuzzy-fmt/figs/'
bindir = '.'
if os.path.isdir('figs'):
    figsdir = 'figs/'
    bindir = '../..'

# always remember to build the executable before running it
#system('scons -U')

def run_homogeneous(n_reduced, temperature, pot = ""):
    nspheres = round(n_reduced*2**(-5.0/2.0)*30**3)
    #width = (nspheres/density)**(1.0/3) # to get density just right!
    filename = '%s/mc%s-%.4f-%.4f' % (figsdir, pot, n_reduced, temperature)
    system("srun --mem=60 -J soft-%.4f-%.4f time nice -19 %s/soft-monte-carlo %d 0.01 0.001 %s.dat periodxyz 30 kT %g potential '%s' > %s.out 2>&1 &" %
           (n_reduced, temperature, bindir, nspheres, filename, temperature, pot, filename))

def run_walls(n_reduced, nspheres, temperature):
    if nspheres == 0:
        nspheres = round(n_reduced*2**(-5.0/2.0)*30**2*32) # reasonable guess
    filename = '%s/mcwalls-%.4f-%.4f-%d' % (figsdir, n_reduced, temperature, nspheres)
    system("srun --mem=60 -J softwalls-%.4f-%.4f time nice -19 %s/soft-monte-carlo %d 0.01 0.001 %s.dat periodxy 30 wallz 30 kT %g > %s.out 2>&1 &" %
           (n_reduced, temperature, bindir, nspheres, filename, temperature, filename))

def run_test_particle(n_reduced, temperature, testp_sigma,testp_eps,pot = ""):
    nspheres = round(n_reduced*2**(-5.0/2.0)*30**3)
    #width = (nspheres/density)**(1.0/3) # to get density just right!
    filename = '%s/mc_testp_%s-%.4f-%.4f' % (figsdir, pot, n_reduced, temperature)
    system("srun --mem=60 -J soft-testp-%.4f-%.4f time nice -19 %s/soft-monte-carlo %d 0.01 0.001 %s.dat periodxyz 30 kT %g TestP %f testp_eps %f potential '%s' > %s.out 2>&1 &" %
           (n_reduced, temperature, bindir, nspheres, filename, temperature, testp_sigma, testp_eps, pot, filename))

#for reduced_density in [0.2, 0.5, 1.0, 1.5]:
#  for reduced_temp in [0.2, 0.5, 1.0]:
#      print 'hello %g %g' % (reduced_density, reduced_temp)
#      run_walls(reduced_density, 0, reduced_temp)


argon_sigma=3.405
argon_eps=119.8

run_homogeneous(.42,.42,"wca")
run_test_particle(0.83890,0.71,argon_sigma,argon_eps,"wca")
run_test_particle(0.957,2.48,argon_sigma,argon_eps,"wca")
run_test_particle(0.957,1.235,argon_sigma,argon_eps,"wca")
run_test_particle(0.5844,1.235,argon_sigma,argon_eps,"wca")
run_test_particle(1.095,2.48,argon_sigma,argon_eps,"wca")




#for reduced_density in [0.7, 0.8, 0.9]:
#    for temp in [10]:
#        run_homogeneous(reduced_density, temp, "wca")
#
#run_homogeneous(0.83890, 0.71) #fig 11
#run_homogeneous(0.957,2.48)    #fig 12
#run_homogeneous(01.095, 2.48)  #fig 13
#run_homogeneous(.5844,1.235)   #fig 14
# run_homogeneous(0.5, 0.0001)
# run_homogeneous(0.6, 0.0001)
# run_homogeneous(0.7, 0.0001)
# run_homogeneous(0.8, 0.0001)

# run_homogeneous(0.1, 0.001)
# run_homogeneous(0.2, 0.001)
# run_homogeneous(0.3, 0.001)
# run_homogeneous(0.4, 0.001)
# run_homogeneous(0.5, 0.001)
# run_homogeneous(0.6, 0.001)
# run_homogeneous(0.7, 0.001)
# run_homogeneous(0.8, 0.001)

# run_homogeneous(0.1, 0.01)
# run_homogeneous(0.2, 0.01)
# run_homogeneous(0.3, 0.01)
# run_homogeneous(0.4, 0.01)
# run_homogeneous(0.5, 0.01)
# run_homogeneous(0.6, 0.01)
# run_homogeneous(0.7, 0.01)
# run_homogeneous(0.8, 0.01)

# run_homogeneous(0.1, 0.1)
# run_homogeneous(0.2, 0.1)
# run_homogeneous(0.3, 0.1)
# run_homogeneous(0.4, 0.1)
# run_homogeneous(0.5, 0.1)
# run_homogeneous(0.6, 0.1)
# run_homogeneous(0.7, 0.1)
# run_homogeneous(0.8, 0.1)

# for reduced_density in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]:
#     for temp in [0.001, 0.01, 0.1, 1.0]:
#         run_homogeneous(reduced_density, temp, "wca")

