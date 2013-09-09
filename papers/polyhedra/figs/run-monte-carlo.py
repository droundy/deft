#!/usr/bin/python
from __future__ import division
from numpy import *
import os

iterations = 100000000
R = sqrt(3)/2
neighborR = 0.5
de_density = 0.01
dr = 0.01

scale = .05
theta_scale = .05


figsdir = 'papers/polyhedra/figs'
bindir = '.'
if os.path.isdir('figs'):
    figsdir = 'figs'
    bindir = '../..'

os.system("scons %s/polyhedra-monte-carlo" %(bindir))

def run_walls(ff, N, dim, shape):
  memory = N/100 # fixme: better guess
  name = "polyhedraMC-walls-%4.2f-%i-%s" %(ff, N, shape)
  filename = 'walls-%4.2f' % ff
  scriptname = "%s/%s-%i-%s.tmp.sh" %(figsdir, filename, N, shape)
  outname = "%s/%s-%i-%s.out" %(bindir, filename, N, shape)
  directory = '%s/mc' % figsdir
  parameters = ["--dir %s" %directory,
                "--filename %s" %filename,
                "--N %i" %N,
                "--iterations %li" %iterations,
                "--periodx %g" %dim,
                "--periody %g" %dim,
                "--wallz %g" %dim,
                "--shape %s" %shape,
                "--R %g" %R,
                "--neighborR %g" %neighborR,
                "--dr %g" %dr,
                "--de_density %g" %de_density,
                "--scale %g" %scale,
                "--theta_scale %g" %theta_scale ]
  command = "time nice -19 %s/polyhedra-monte-carlo" %bindir
  script = open(scriptname, 'w')
  script.write("#!/bin/bash\n")
  script.write("#SBATCH --mem-per-cpu=%i\n" % memory)
  script.write("##SBATCH --mail-type ALL\n")
  script.write("##SBATCH --mail-user paho@paholg.com\n")
  script.write("#SBATCH --output %s\n\n" % outname)

  script.write("echo \"Starting polyhedra-monte-carlo with estimated memory use: %i.\"\n\n" %memory)
  script.write(command)
  for param in parameters:
    script.write(" " + param)
  script.write("\n")

  script.close()

  os.system("sbatch -J %s %s\n" %(name, scriptname))



run_walls(.3, 2400, 20, 'cube')
run_walls(.4, 3200, 20, 'cube')
run_walls(.5, 4000, 20, 'cube')
run_walls(.6, 4800, 20, 'cube')
run_walls(.7, 5600, 20, 'cube')
run_walls(.73, 5832, 20, 'cube')
run_walls(.30, 2197, 20, 'truncated_tetrahedron')
run_walls(.37, 2744, 20, 'truncated_tetrahedron')
run_walls(.46, 3375, 20, 'truncated_tetrahedron')
