#!/usr/bin/python2

import re, string, sys, os

import facfile

if len(sys.argv) > 1:
    texfs = sys.argv[1:]
else:
    print("Need arguments")
    os.exit(1)

documentclassre = re.compile(r'\\documentclass(\[[^\]]+\])?{')
graphicre = re.compile(r'^\s*\\includegraphics(\[[^\]]+\])?{([^}]+)}', re.MULTILINE)
inputre = re.compile(r'\\input(\[[^\]]+\])?{([^}]+)}')

for t in texfs:
    fac = facfile.facfile(os.path.dirname(t)+'/.tex.fac')

    fname = os.path.basename(t)
    f = open(t, 'r')
    latex = f.read()
    inputs = {fname}
    outputs = {fname[:-3]+'pdf'}

    os.chdir(os.path.dirname(t)) # so we can look up figures easily
    def graphics_name(x):
        if len(x) > 4 and x[-4:] in ['.pdf', '.png', '.jpg']:
            return x
        return x+'.pdf'
    inputs |= set([graphics_name(x[1]) for x in graphicre.findall(latex)])
    inputs |= set([x[1]+'.tex' for x in inputre.findall(latex)])
    fac.default('pdflatex %s && pdflatex %s && bibtex %s && pdflatex %s'
                % (fname, fname, fname[:-4], fname), inputs, outputs)


