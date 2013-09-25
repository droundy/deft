#!/usr/bin/python
from __future__ import division
from numpy import *
import os, sys, argparse

import arguments

iterations = 100000000
neighborR = 0.5
de_density = 0.01
dr = 0.01
dim = 20

N = 1000

scale = .05
theta_scale = .05


figsdir = 'papers/polyhedra/figs'
bindir = '.'
if os.path.isdir('figs'):
    figsdir = 'figs'
    bindir = '../..'

os.system("scons %s/polyhedra-monte-carlo" %(bindir))

def run_mc(ff, shape, celltype="periodic", ratio=1):
  memory = N/100 # fixme: better guess
  name = "polyhedraMC-%s-%4.2f-%i-%s" %(celltype, ff, N, shape)
  filename = '%s-%4.2f' %(celltype, ff)
  scriptname = "%s/%s-%i-%s.tmp.sh" %(figsdir, filename, N, shape)
  outname = "%s/%s-%i-%s.out" %(bindir, filename, N, shape)
  directory = '%s/mc' % figsdir
  if celltype == "periodic":
    cellparam = "--periodz %g" %dim
  elif celltype == "walls":
    cellparam = "--wallz %g" %dim
  else:
    print("invalid cell type")
    exit(1)
  parameters = [
    "--dir %s" %directory,
    "--filename %s" %filename,
    "--N %i" %N,
    "--iterations %li" %iterations,
    "--periodx %g" %dim,
    "--periody %g" %dim,
    cellparam,
    "--shape %s" %shape,
    "--ff %g" %ff,
    "--neighborR %g" %neighborR,
    "--dr %g" %dr,
    "--de_density %g" %de_density,
    "--ratio %g" %ratio
  ]

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

if(len(sys.argv) == 1):
  run_mc(.3, "cube")
  run_mc(.5, "cube")
  run_mc(.7, "cube")
  run_mc(.3, "cuboid", ratio=2)
  run_mc(.5, "cuboid", ratio=2)
  run_mc(.7, "cuboid", ratio=2)

else:
  parser = argparse.ArgumentParser(
    description='Plot density of polyhedra.', parents = [arguments.parser])

  args = parser.parse_args()

  ff = args.ff
  polyhedron = args.shape
  ratio = args.ratio

  if args.periodic:
    celltype = 'periodic'
  else:
    celltype = 'walls'

  run_mc(ff, polyhedron, ratio=ratio, celltype=celltype)
