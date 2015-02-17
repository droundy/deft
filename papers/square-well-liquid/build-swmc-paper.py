#!/usr/bin/env python3

import os, sys
import subprocess as sp

# simulation parameters
ww = 1.3
ff = 0.3
N = 20
iterations = 1e7

# directory names, etc.
swmc_dir = os.path.dirname(os.path.realpath(__file__))
fig_dir = os.path.join(swmc_dir,"figs")
swmc = "square-well-monte-carlo"
script_env = "python2"

def fig_name(fig):
    return os.path.join(fig_dir,'plot-'+fig+'.py')

# get stuff inside brackets next to scons hook for an argement
def get_field(script_name,field):
    with open(script_name) as script:
        for line in script:
            if "#arg %s" % field in line:
                return ''.join(line.split()[3:])[1:-1]

# execute processes from a list of arguments, adding them to a process list
ps = []
def execute(args):
    print(" ".join(args))
    ps.append(sp.Popen(args))

# generate all figures
os.chdir(swmc_dir)
for fig in ["histograms","dos","samples","uhc","scaling","transitions"]:

    script_name = fig_name(fig)

    if fig != "scaling":
        Ns = str(N)
    else:
        Ns = get_field(script_name,"all_Ns")

    if fig != "transitions":
        versions = [get_field(script_name,"versions")]
    else:
        versions = [get_field(script_name,"method1").replace("'",""),
                    get_field(script_name,"method2").replace("'","")]

    execute([script_env,script_name,str(ww),str(ff),Ns] + versions)

# wait for all figure scripts to finish
for p in ps:
    p.wait()

# build histogram-paper.tex
execute(["pdflatex","-interaction=nonstopmode","-recorder","histogram-paper.tex"])
