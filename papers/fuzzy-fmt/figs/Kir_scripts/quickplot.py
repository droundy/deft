#!/usr/bin/python2

import os
import numpy as np
import matplotlib.pyplot as plt

#for gwidth=0.001
os.system('figs/new-melting.mkdat --kT 2 --n 1.08 --d carrots_kT2_n1.08_gw0.001 --fv -1 --gw .001  --dx 0.01')
os.system('cd carrots_kT2_n1.08_gw0.001')
os.system('cat *alldat.dat >> plot.dat')
os.system('cd ..')
os.system('./new-melting_anyplot_script.py carrots_kT2_n1.08_gw0.001 --ftemp 2 --ydiff --xfv')

#for gwidth=0.01
#os.system('figs/new-melting.mkdat --kT 2 --n 1.08 --d carrots_kT2_n1.08_gw0.01 --fv -1 --gw .01  --dx 0.01')
#os.system('cd carrots_kT2_n1.08_gw0.01')
#os.system('cat *alldat.dat >> plot.dat')
#os.system('cd ..')
#os.system('./new-melting_anyplot_script.py carrots_kT2_n1.08_gw0.01 --ftemp 2 --ydiff --xfv')

#for gwidth=0.1
#os.system('figs/new-melting.mkdat --kT 2 --n 1.08 --d carrots_kT2_n1.08_gw0.1 --fv -1 --gw .1  --dx 0.01')
#os.system('cd carrots_kT2_n1.08_gw0.1')
#os.system('cat *alldat.dat >> plot.dat')
#os.system('cd ..')
#os.system('./new-melting_anyplot_script.py carrots_kT2_n1.08_gw0.1 --ftemp 2 --ydiff --xfv')

#for gwidth=0.4
#os.system('figs/new-melting.mkdat --kT 2 --n 1.08 --d carrots_kT2_n1.08_gw0.4 --fv -1 --gw .4  --dx 0.01')
#os.system('cd carrots_kT2_n1.08_gw0.4')
#os.system('cat *alldat.dat >> plot.dat')
#os.system('cd ..')
#os.system('./new-melting_anyplot_script.py carrots_kT2_n1.08_gw0.4 --ftemp 2 --ydiff --xfv')
