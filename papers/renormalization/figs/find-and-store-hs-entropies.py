#!/usr/bin/python2
from __future__ import division
import numpy as np
import gatherandcalculate, os

ww = 1.3
L = 2.84

i_values = [1,2]
for i in i_values:
    dbase = 'data/scrunched-ww%04.2f-L%04.2f/i%01d' % (ww, L,  i)
    os.system('rm -f %s/*/*/Sexc.dat' % dbase)
    S, Ns = gatherandcalculate.Sexc_hardsphere_Ns(dbase)
    for n in range(len(S)):
        # here we save this single numb er for later...
        np.savetxt(dbase+'/N%03d/absolute/Sexc.dat' % Ns[n], [S[n]], fmt='%.9g')
