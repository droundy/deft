#!/usr/bin/python2
#Run this program from /deft/papers/fuzzy-fmt by entering ./nm_meanFE_vs_gw_data.py [filename.dat]
#Look for a filename.dat (file with data to plot) with a name like:  
#FE_vs_gw_kT2_n1.3_fv0_dx0.5_mcerror0.001_mcconstant5_mcprefactor50000_seeds10.dat


import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Creates a data file with meanFE and uncertainty.')
parser.add_argument('datafile', metavar='datafile', type=str,
                    help='file with data to process') 
args=parser.parse_args()

datafile= args.datafile      

os.system('./nm_plot_FE_vs_gw.py %s >> mean%s' % (datafile, datafile))
