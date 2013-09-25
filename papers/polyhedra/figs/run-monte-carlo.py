#!/usr/bin/python
from __future__ import division
from numpy import *
import os, sys, argparse

import arguments

iterations = 10000000
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
  memory = N/40 # fixme: better guess
  filename = '%s-%4.2f' %(celltype, ff)
  if ratio == 1: ratio_name = ""
  else: ratio_name = "_%05.2f" %ratio
  name = "poly-%s-%s%s-%i" %(filename, shape, ratio_name, N)
  scriptname = "%s/%s-%s%s-%i.tmp.sh" %(figsdir, filename, shape, ratio_name, N)
  outname = "%s/%s-%s%s-%i.out" %(bindir, filename, shape, ratio_name, N)
  if celltype == "periodic":
    cellparam = "--periodz %g" %(dim*ratio)
  elif celltype == "walls":
    cellparam = "--wallz %g" %(dim*ratio)
  else:
    print("invalid cell type")
    exit(1)
  parameters = [
    "--N %i" %N,
    "--iterations %li" %iterations,
    "--periodx %g" %dim,
    "--periody %g" %dim,
    cellparam,
    "--shape %s" %shape,
    "--ff %g" %ff,
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
  run_mc(.4, "cube")
  run_mc(.5, "cube")
  run_mc(.6, "cube")
  run_mc(.7, "cube")
  run_mc(.3, "cuboid", ratio=2)
  run_mc(.4, "cuboid", ratio=2)
  run_mc(.5, "cuboid", ratio=2)
  run_mc(.6, "cuboid", ratio=2)
  run_mc(.7, "cuboid", ratio=2)

else:
  parser = argparse.ArgumentParser(
    description='Create batch script and run monte carlo simulation(s) using sbatch.', parents = [arguments.parser])

  args = parser.parse_args()

  ff = args.ff
  polyhedron = args.shape
  if args.ratio == 0: ratio = 1
  else: ratio = args.ratio

  if args.periodic:
    celltype = 'periodic'
  else:
    celltype = 'walls'

  run_mc(ff, polyhedron, ratio=ratio, celltype=celltype)
