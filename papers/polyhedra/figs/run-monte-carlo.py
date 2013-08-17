#!/usr/bin/python
from __future__ import division
import os

iterations = 100000000
R = 1.7320508075688772
dim = 20

scale = .005
theta_scale = .05


figsdir = 'papers/polyhedra/figs'
bindir = '.'
if os.path.isdir('figs'):
    figsdir = 'figs'
    bindir = '../..'

os.system("time scons %s/polyhedra-monte-carlo" %(bindir))

memory = 200 # fixme: better guess
def run_walls(ff, N, shape):
  scriptname = "%s/polyhedraMC-walls-%4.2f-%i-%s.tmp.sh" %(figsdir, ff, N, shape)
  outname = "%s/polyhedraMC-walls-%4.2f-%i-%s.out" %(bindir, ff, N, shape)
  filename = '%s/mc/polyhedraMC-walls-%4.2f' % (figsdir, ff)
  jobID = "polyhedraMC-walls-%4.2f-%i-%s" %(ff, N, shape)
  command = "time nice -19 %s/polyhedra-monte-carlo %i %i %s wallz periodx periody \
R %g dimensions %g %g %g scale %g theta_scale %g" %(bindir, N, iterations, filename,
                                                    R, dim, dim, dim, scale, theta_scale)
  script = open(scriptname, 'w')
  script.write("#SBATCH --mem-per-cpu=%i\n" %memory)
  script.write("##SBATCH --mail-type ALL\n")
  script.write("##SBATCH --mail-user paho@paholg.com\n")
  script.write("#SBATCH --output %s\n\n" %outname)

  script.write("echo \"Starting polyhedra-monte-carlo with estimated memory use: %i.\"\n\n" %memory)
  script.write("%s\n" %(command))

  os.system("sbatch -J %s %s\n" %(jobID, scriptname))
  os.system("rm %s" %(scriptname))




run_walls(.6, 600, 'cube')
