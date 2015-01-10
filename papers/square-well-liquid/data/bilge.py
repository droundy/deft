#!/usr/bin/python2

for (ww, ff, N, method) in [(1.3, 0.3, 20, 'nw'),
                            (1.3, 0.3, 20, 'flat'),
                            (1.3, 0.3, 20, 'gaussian'),
                            (1.3, 0.3, 20, 'wang_landau'),
                            (1.3, 0.3, 20, 'walkers'),
                            (1.3, 0.3, 20, 'kT 1'),
                            (1.3, 0.3, 20, 'kT 2')]:
    outputs = ["periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat"
               % (ww, ff, N, method.replace(' ',''), postfix)
               for postfix in ['g', 'E', 'lnw', 'rt']]

    print '| cd ../../.. && ./square-well-monte-carlo --%s --N %d --initialize 1000 --ff %g --ww %g --iterations 10000' % (method, N, ff, ww)
    print '< ../../../square-well-monte-carlo'

    for o in outputs:
        print '>', o
    print
