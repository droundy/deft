#!/usr/bin/python2

import glob

for py in glob.glob('*/*/*/*.py') + glob.glob('*/*/*.py') + glob.glob('*/*.py') + glob.glob('*.py'):
    print '| python -m compileall', py
    print '>', py+'c'
    print

