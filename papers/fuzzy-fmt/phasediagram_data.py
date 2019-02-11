#!/usr/bin/python2
#NOTE: Run this plot script from directory deft/papers/fuzzy-fmt 
#with comand ./phasediagram_data.py --tar [optional]

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os, glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--tar', action='store_true',
                    help='--tar for use tarazona tensor')
args=parser.parse_args()

for kT in np.arange(0.05, 1.25, 0.05):
    #print (kT)
    #cmd = 'rq run -J plot-pressure.py --kT %g >> phasenew.dat' % (kT)
    if args.tar :
        cmd = './plot-pressure.py --kT %g --tar >> phasenewtar.dat' % (kT)
    else :
        cmd = './plot-pressure.py --kT %g >> phasenew.dat' % (kT)
    print(cmd)
    os.system(cmd)
