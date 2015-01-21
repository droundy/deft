#!/usr/bin/python2

for method in ["nw","bubble_suppression","robustly_optimistic","gaussian","wang_landau","walker_optimization"]:
    for ww in [1.3]:
        for ff in [0.3]:
            for N in [5,6,7,8,9,10,20]:
                outputs = ["periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat"
                           % (ww, ff, N, method.replace(' ',''), postfix)
                           for postfix in ['g', 'E', 'lnw', 's']]

                print '| cd ../../.. && ./square-well-monte-carlo --%s --N %d --ff %g --ww %g --iterations 1000000' % (method, N, ff, ww)
                print '< ../../../square-well-monte-carlo'

                for o in outputs:
                    print '>', o
                print
