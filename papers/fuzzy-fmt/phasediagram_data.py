#!/usr/bin/python2

#This program generates the phasenew.dat (or phasenewtensor.dat) file 
#used by plot-phasediagram to plot phase diagrams. It runs
#plot-pressure.py for various temperatures and directs the 
#printed output into the phasenew.dat (or phasenewtensor.dat) file.

#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./phasediagram_data.py [OPTIONAL: --tensor]

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tensor', action='store_true',
                    help='use tensor weight')
args=parser.parse_args()

for kT in np.arange(0.05, 1.25, 0.05):
    if args.tensor :
        cmd = 'rq run -J plot-pressure.py --kT %g --tensor >> phasenewtensor.dat' % (kT)
    else :
        cmd = 'rq run -J plot-pressure.py --kT %g >> phasenew.dat' % (kT)
    print(cmd)
    os.system(cmd)
