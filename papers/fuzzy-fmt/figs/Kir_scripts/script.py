#!/usr/bin/python2

import os

os.system('rm newmeltdataout.dat')

for gwidth in [.46, .459, .458, .4578, .4576, .4574, .4572, .457, .456]:
    os.system('../new-melting.mkdat 1.3 .1 %g 2' % (gwidth))

os.system('g++ -o newmeltdataout newmeltdataout.cpp')
os.system('./newmeltdataout')
os.system('gnuplot newmeltdataout.gnu')
