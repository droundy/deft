#!/usr/bin/python3

import facfile
import string

src = facfile.facfile('.monte-carlos.fac')

monte_carlos = """
   monte-carlo soft-monte-carlo pair-monte-carlo
   triplet-monte-carlo polyhedra-monte-carlo polyhedra-talk
   square-well-monte-carlo
   radial-distribution-monte-carlo
   free-energy-monte-carlo free-energy-monte-carlo-infinite-case
""".split()

utility_files = """
   polyhedra square-well
""".split()

for x in utility_files:
    src.compile('src/Monte-Carlo/%s.cpp' % x)

for x in monte_carlos:
    src.link('src/Monte-Carlo/%s.cpp' % x, x)

def add_parameters(method):
    if method in ['tmmc', 'oetmmc']:
        return method + ' --min_samples 10000'
    return method

T_sims = ["kT %g" %kT for kT in [i*.1 for i in range(1,10)] + range(1,10)]
hist_methods = ['simple_flat','wang_landau','tmmc','oetmmc']
for i in range(len(hist_methods)):
    hist_methods.append(hist_methods[i]+'_oe')

for method in ["nw"] + T_sims + hist_methods:
    for ww in [1.3, 1.5]:
        for ff in [0.3, 0.8]:
            for N in range(5,31):
                outputs = ["papers/histogram/data/periodic-ww%04.2f-ff%04.2f-N%i-%s-%s.dat"
                           % (ww, ff, N, method.replace(' ',''), postfix)
                           for postfix in ['E','lnw','transitions','os','ps','g']],
                src.rule('./square-well-monte-carlo --N %d --%s --ff %g --ww %g --iterations 3000000'
                         % (N, add_parameters(method.replace('_oe',' --optimized_ensemble')), ff, ww),
                         ['square-well-monte-carlo'], outputs)
