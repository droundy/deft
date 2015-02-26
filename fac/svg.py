#!/usr/bin/python3

import os, glob

os.chdir('..')

for svg in glob.glob('papers/*/figs/*.svg') + glob.glob('papers/*/*.svg'):
    print("? inkscape --export-eps %s %s" % (svg[:-3]+'eps', svg))
    print(">", svg[:-3]+'eps')
    print("C", os.getenv("HOME")+'/.config/inkscape')
    print()
    print("? epstopdf --outfile %s %s" % (svg[:-3]+'pdf', svg[:-3]+'eps'))
    print(">", svg[:-3]+'pdf')
    print("<", svg[:-3]+'eps')
    print()
