#!/usr/bin/python2
#NOTE: Run this script from deft/papers/fuzzy-fmt with the 
#command ./nm_hist_displayplots.py

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import os

for mcprefactor in [10000, 20000, 30000, 40000, 50000]: 
  os.system('./nm_plot_histdata.py Histogram_%g_mcerror0.001/Hist*%g_gw*.dat &' % (mcprefactor, mcprefactor)) 
  #the & symbol runs it in the background so all the plots are displayed without pauses due to plt.show() command in nm_plot_histdata.py

