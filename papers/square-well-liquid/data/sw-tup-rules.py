#!/usr/bin/python

import sys

name, ww, ff, N, method = sys.argv
ww = float(ww)
ff = float(ff)
N = int(N)

outputs = ["periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat"
           % (ww, ff, N, method.replace(' ',''), postfix)
           for postfix in ['g', 'E', 'lnw', 'rt']]

cmd = 'cd ../../.. && ./square-well-monte-carlo --%s --N %d --initialize 1000 --ff %g --ww %g --iterations 10000' % (method, N, ff, ww)

rule = ': ../../../square-well-monte-carlo |> %s |>' % cmd

for o in outputs:
    rule += ' ' + o

print rule
