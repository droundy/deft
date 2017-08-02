#!/usr/bin/env python2

import os, sys

ww = float(sys.argv[1])
ff = float(sys.argv[2])
lenx = float(sys.argv[3])
lenyz = float(sys.argv[4])
min_T = float(sys.argv[5])
method_name = sys.argv[6]

fname = 'ww%.2f-ff%.2f-%gx%g-%s' % (ww, ff, lenx, lenyz, method_name)

os.system('sbatch -J lv:%s -o data/lv/%s.out data/run-sticky-wall.py '
           % (fname,fname) + ' '.join(sys.argv[1:]))
