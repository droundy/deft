#!/usr/bin/python
from __future__ import division
import os


iterations = 100000000
start_keeping = 10000
R = 1.7320508075688772
dim = 20

scale = .005
theta_scale = .05


figsdir = 'papers/polyhedra/figs'
bindir = '.'
if os.path.isdir('figs'):
    figsdir = 'figs/'
    bindir = '../..'

os.system("time scons polyhedra-monte-carlo")

memory = 20 # fixme: better guess
def run_walls(ff, N, shape):
    fname = '%s/mc/wallsMC-%4.2f' % (figsdir, ff)
    jobID = "polyhedra-%4.2f-%i-%s" %(ff, N, shape)
    command = "time nice -19 %s/polyhedra-monte-carlo %i %i %s wallz periodx periody R %g dimensions %g %g %g start_keeping %i scale %g theta_scale %g" %(bindir, N, iterations, fname, R, dim, dim, dim, start_keeping, scale, theta_scale)
    os.system("srun --mem=%g -J %s %s" %(memory, jobID, command))




run_walls(.6, 600, 'cube')
