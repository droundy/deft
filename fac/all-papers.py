#!/usr/bin/python3

import os

papers = """
   histogram pair-correlation contact hard-sphere-free-energy
   square-well-fluid
   fuzzy-fmt water-saft polyhedra electrostatics renormalization
""".split()

for paper in papers:
  print("| python3 latex.py papers/"+paper+"/paper.tex")
  print("< config.py")
  print("> ../papers/"+paper+"/.paper.tex.fac")
  print("c .pyc\n")

  print("| cd ../papers/%s && python2 ../../fac/figs.py figs/*.py > .figs.fac"
        % paper)
  print("> ../papers/"+paper+"/.figs.fac")
  print("c .pyc\n")

for project in ['thesis-roth']:
  print("| python3 latex.py papers/"+project+"/project.tex")
  print("< config.py")
  print("> ../papers/"+project+"/.project.tex.fac")
  print("c .pyc\n")

  print("| cd ../papers/%s && python2 ../../fac/figs.py figs/*.py > .figs.fac"
        % project)
  print("> ../papers/"+project+"/.figs.fac")
  print("c .pyc\n")

# Add thesis-vischer/Thesis.tex here when figures are available
other_texs = """
   histogram/thesis/thesis.tex
""".split()

for tex in other_texs:
  d = os.path.dirname(tex)
  facd = os.getcwd()
  print("| python3 latex.py papers/"+tex)
  print("< config.py")
  print("> ../papers/"+d+"/.%s.fac" % os.path.basename(tex))
  print("c .pyc\n")

  print("| cd ../papers/%s && python2 %s/figs.py figs/*.py > .figs.fac"
        % (d, facd))
  print("> ../papers/"+d+"/.figs.fac")
  print("c .pyc\n")
