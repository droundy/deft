#!/usr/bin/python2

import glob, re, string, sys

if len(sys.argv) > 1:
    texfs = sys.argv[1:]
else:
    texfs = glob.glob('*.tex')

documentclassre = re.compile(r'\\documentclass(\[[^\]]+\])?{')
graphicre = re.compile(r'^\s*\\includegraphics(\[[^\]]+\])?{([^}]+)}', re.MULTILINE)
inputre = re.compile(r'\\input(\[[^\]]+\])?{([^}]+)}')

for fname in texfs:
    f = open(fname, 'r')
    latex = f.read()
    print '| pdflatex %s && pdflatex %s && bibtex %s && pdflatex %s' % (fname, fname, fname[:-4], fname)
    print '>', fname[:-3]+'bbl'
    print '>', fname[:-3]+'pdf'
    if documentclassre.search(latex):
        for pdf in [x[1]+'.pdf' for x in graphicre.findall(latex)]:
            print '<', pdf
        for tex in [x[1]+'.tex' for x in inputre.findall(latex)]:
            print '<', tex
    print

