#!/usr/bin/python2
from __future__ import division
import matplotlib, sys
if 'show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy, time, os

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import readandcompute

if len(sys.argv) != 6:
    print "usage:  python2 %s ww ff N minT iterations_per_frame" % sys.argv[0]
    exit(1)

if os.path.exists('paper.tex'):
    os.chdir('../..')

assert os.path.exists('SConstruct')
# switch to deft project directory and build SWMC
assert not os.system('fac square-well-monte-carlo')

ww = float(sys.argv[1])
ff = float(sys.argv[2])
N = int(sys.argv[3])
min_T = float(sys.argv[4])
iterations_per_frame = int(float(sys.argv[5]))

mem_estimate = 10 + 0.1*N # it actually also depends on ww, but I'm ignoring that for now.

datadir = 'papers/square-well-fluid/data/mc'
fname = 'ww%.2f-ff%.2f-N%d' % (ww, ff, N)
moviedir = datadir + '/movies/ww%.2f-ff%.2f-N%d' % (ww, ff, N)

os.system('mkdir -p ' + datadir)
os.system('rm -rf ' + moviedir)
os.system('mkdir -p ' + moviedir)

if os.system('which srun'):
    # srun doesn't exist so just run on this machine
    cmd = "time nice -19 ./square-well-monte-carlo"
else:
    cmd = "srun --mem=%d -J %s time nice -19 ./square-well-monte-carlo" % (mem_estimate, fname)

cmd += " --ww %g --ff %g --N %d" % (ww, ff, N)

cmd += " --iterations 1 --init_iters %d --golden" % (iterations_per_frame)

cmd += ' --de_g 0.01' # nice high-resolution radial distribution function data

cmd += ' --min_T %g' % min_T

cmd += ' --data_dir %s --filename %s' % (datadir, fname)

fig, ax = plt.subplots()

total_iterations = 0
for frame in xrange(1000000):
    assert not os.system(cmd)
    total_iterations += iterations_per_frame

    plt.cla()

    basename = datadir+'/ww%.2f-ff%.2f-N%d' % (ww,ff,N)
    e, diff = readandcompute.e_diffusion_estimate(basename)

    T, u, cv, s, minT = readandcompute.T_u_cv_s_minT(basename)
    ax.plot(u/N, T, 'k-')
    ax.set_ylim(0, 3)
    ax.axhline(minT, color='r', linestyle=':')

    ax.axvline(-readandcompute.max_entropy_state(basename)/N, color='r', linestyle=':')
    ax.axvline(-readandcompute.min_important_energy(basename)/N, color='r', linestyle='--')

    e, init_hist = readandcompute.e_and_total_init_histogram(basename)
    ax.plot(e, init_hist, 'k-')


    ax.set_xlabel(r'$E/N$')
    ax.set_ylim(0, init_hist.max()*1.1)
    ax.set_ylabel(r'$N$')
    plt.title(r'Initialization histogram with $\lambda = %g$, $\eta = %g$, and $N = %d$ (%.0fM iters)'
              % (ww, ff, N, total_iterations/1.0e6))

    fname = moviedir+'/frame%06d.png' % (frame)
    print 'saving', fname
    plt.savefig(fname)
    os.system("convert -delay 100 %s/frame*.png %s/movie.gif"
              % (moviedir, moviedir)) # make the movie

