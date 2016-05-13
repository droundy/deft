import scipy as sp
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import pylab as plt
import matplotlib
import RG
import SW
import numpy as np
import time
import integrate
import os
import sys

temp = plt.linspace(0.6,1.28,20)

for T in temp:
        os.system('python RG_plots_test.py %f'%T)


#                                                                                             #
#                                                                                             #
######################################## END PLOTTING STUFF ###################################
