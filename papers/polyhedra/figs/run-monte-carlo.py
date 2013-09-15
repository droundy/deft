#!/usr/bin/python
from __future__ import division
from numpy import *
import os

iterations = 100000000
R = 1
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

def run_mc(ff, N, dim, shape, celltype="periodic"):
  memory = N/100 # fixme: better guess
  name = "polyhedraMC-%s-%4.2f-%i-%s" %(celltype, ff, N, shape)
  filename = '%s-%4.2f' %(celltype, ff)
  scriptname = "%s/%s-%i-%s.tmp.sh" %(figsdir, filename, N, shape)
  outname = "%s/%s-%i-%s.out" %(bindir, filename, N, shape)
  directory = '%s/mc' % figsdir
  if celltype == "periodic":
    cellparam = "--periodz %g" %dim
  else:
    cellparam = "--wallz %g" %dim
  parameters = ["--dir %s" %directory,
                "--filename %s" %filename,
                "--N %i" %N,
                "--iterations %li" %iterations,
                "--periodx %g" %dim,
                "--periody %g" %dim,
                cellparam,
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


run_mc(.58, 2744, 20, 'truncated_tetrahedron', 'periodic')
run_mc(.46, 2197, 20, 'truncated_tetrahedron', 'periodic')
run_mc(.36, 1728, 20, 'truncated_tetrahedron', 'periodic')

