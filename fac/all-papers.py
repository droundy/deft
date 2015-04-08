#!/usr/bin/python3

import os

papers = """
   histogram pair-correlation contact hard-sphere-free-energy
   fuzzy-fmt water-saft polyhedra electrostatics renormalization
""".split()

for paper in papers:
  print("| python3 latex.py papers/"+paper+"/paper.tex")
  print("< config.py")
  print("> ../papers/"+paper+"/.tex.fac")
  print("c .pyc\n")

  print("| cd ../papers/%s && python2 ../../fac/figs.py figs/*.py > .figs.fac"
        % paper)
  print("> ../papers/"+paper+"/.figs.fac")
  print("c .pyc\n")
