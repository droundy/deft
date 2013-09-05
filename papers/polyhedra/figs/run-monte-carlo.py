#!/usr/bin/python
from __future__ import division
from numpy import *
import os

iterations = 100000000
R = sqrt(3)/2
dim = 20
de_density = 0.01

scale = .05
theta_scale = .05


figsdir = 'papers/polyhedra/figs'
bindir = '.'
if os.path.isdir('figs'):
    figsdir = 'figs'
    bindir = '../..'

os.system("scons %s/polyhedra-monte-carlo" %(bindir))

def run_walls(ff, N, shape):
  memory = N/100 # fixme: better guess
  name = "polyhedraMC-walls-%4.2f-%i-%s" %(ff, N, shape)
  filename = 'walls-%4.2f' % ff
  scriptname = "%s/%s-%i-%s.tmp.sh" %(figsdir, filename, N, shape)
  outname = "%s/%s.out" %(bindir, name)
  directory = '%s/mc' % figsdir
  command = "time nice -19 %s/polyhedra-monte-carlo %i %i %s %s periodx periody wallz\
 R %g dimensions %g %g %g scale %g theta_scale %g shape %s de_density %g" \
      %(bindir, N, iterations, directory, filename, R, dim, dim, dim, scale, theta_scale,
        shape, de_density)
  script = open(scriptname, 'w')
  script.write("#!/bin/bash\n")
  script.write("#SBATCH --mem-per-cpu=%i\n" % memory)
  script.write("##SBATCH --mail-type ALL\n")
  script.write("##SBATCH --mail-user paho@paholg.com\n")
  script.write("#SBATCH --output %s\n\n" % outname)

  script.write("echo \"Starting polyhedra-monte-carlo with estimated memory use: %i.\"\n\n" %memory)
  script.write("%s\n" %(command))

  script.close()

  os.system("sbatch -J %s %s\n" %(name, scriptname))



run_walls(.3, 2400, 'cube')
run_walls(.4, 3200, 'cube')
run_walls(.5, 4000, 'cube')
run_walls(.6, 4800, 'cube')
run_walls(.7, 5600, 'cube')
run_walls(.8, 6400, 'cube')
