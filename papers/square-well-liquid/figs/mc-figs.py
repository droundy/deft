#!/usr/bin/python2

import sys, string

name, ww, ff, N = sys.argv
ww = float(ww)
ff = float(ff)
N = int(N)

methods = ["nw", "wang_landau", "gaussian", "flat", "walkers", "kT2", "kT1"]

inputs = ["../data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat"
          % (ww, ff, N, method.replace(' ',''), postfix)
          for postfix in ['g', 'E', 'lnw', 'rt']
          for method in methods] + ['styles.pyc']

cmd = 'cd .. && python2 figs/plot-dos.py %f %f %d "%s"' % (ww, ff, N, methods)
outputs = ["periodic-ww%02.0f-ff%02.0f-N%i-dos.pdf" % (ww*100, ff*100, N),
           "periodic-ww%02.0f-ff%02.0f-N%i-weights.pdf" % (ww*100, ff*100, N)]
print ': %s |> %s |> %s' % (string.join(inputs), cmd, string.join(outputs))

cmd = 'cd .. && python2 figs/plot-histograms.py %f %f %d "%s"' % (ww, ff, N, methods)
outputs = ["periodic-ww%02.0f-ff%02.0f-N%i-E.pdf" % (ww*100, ff*100, N)]
print ': %s |> %s |> %s' % (string.join(inputs), cmd, string.join(outputs))

cmd = 'cd .. && python2 figs/plot-uhc.py %f %f %d "%s"' % (ww, ff, N, methods)
outputs = ["periodic-ww%02.0f-ff%02.0f-N%i-u.pdf" % (ww*100, ff*100, N),
           "periodic-ww%02.0f-ff%02.0f-N%i-hc.pdf" % (ww*100, ff*100, N)]
print ': %s |> %s |> %s' % (string.join(inputs), cmd, string.join(outputs))

cmd = 'cd .. && python2 figs/plot-round-trips.py %f %f %d "%s"' % (ww, ff, N, methods)
outputs = ["periodic-ww%02.0f-ff%02.0f-N%i-rt.pdf" % (ww*100, ff*100, N),
           "periodic-ww%02.0f-ff%02.0f-N%i-rt-rate.pdf" % (ww*100, ff*100, N)]
print ': %s |> %s |> %s' % (string.join(inputs), cmd, string.join(outputs))