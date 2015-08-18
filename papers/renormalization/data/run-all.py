#!/usr/bin/python

import sys, os

i = eval(sys.argv[1])
#RG recursion level

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float(eval((sys.argv[3])))
#arg L = [5.0]

Ns = eval(sys.argv[4])
#arg Ns = [range(2,10)]

cmd_fe = 'python run-absolute.py %d %g %d "%s"' % (i,ww,L,Ns)
print(cmd_fe)
os.system(cmd_fe)

cmd_mc = "python run-monte-carlo.py %d %g %d '%s'" %(i,ww,L,Ns)
print(cmd_mc)
os.system(cmd_mc)
