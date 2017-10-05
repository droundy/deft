#!/usr/bin/env python3

import os, glob

if os.system('git diff --exit-code .'):
    print('There are uncommitted changes in this directory.  Commit these before formatting!')
    exit(0)

for f in glob.glob('*.cpp') + glob.glob('*/*.cpp'):
    print('formatting', f)
    os.system('astyle --style=google --indent=spaces=2 %s' % f)
