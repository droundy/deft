#!/usr/bin/python3

import os

papers = """
   histogram pair-correlation contact hard-sphere-free-energy
   square-well-fluid grand
   fuzzy-fmt water-saft polyhedra electrostatics renormalization
""".split()

for paper in papers:
  print("| python3 latex.py papers/"+paper+"/paper.tex")
  print("< config.py")
  print("> ../papers/"+paper+"/.paper.tex.fac")
  print("c .pyc\n")

  print("| cd ../papers/%s && python3 ../../fac/figs.py figs/*.py > .figs.fac"
        % paper)
  print("> ../papers/"+paper+"/.figs.fac")
  print("c .pyc\n")

figs_done = set()

other_texs = """
   thesis-vischer/Thesis.tex
   histogram/sad.tex
   histogram/sad2d.tex
   histogram/ising.tex
   ising/ising-sad.tex
   thesis-kirstie/thesis.tex
""".split()

for tex in other_texs:
  d = os.path.dirname(tex)
  facd = os.getcwd()
  print("| python3 latex.py papers/"+tex)
  print("< config.py")
  print("> ../papers/"+d+"/.%s.fac" % os.path.basename(tex))
  print("c .pyc\n")

  if d not in papers:
    print("* cd ../papers/%s && python3 %s/figs.py figs/*.py > .figs.fac"
          % (d, facd))
    print("> ../papers/"+d+"/.figs.fac")
    print("c .pyc\n")
