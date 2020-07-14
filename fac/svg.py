#!/usr/bin/python3

import os, glob

os.chdir('..')

for svg in glob.glob('papers/*/figs/*.svg') + glob.glob('papers/*/*.svg'):
    print("? cairosvg -o %s %s" % (svg[:-3]+'pdf', svg))
    print(">", svg[:-3]+'pdf')
