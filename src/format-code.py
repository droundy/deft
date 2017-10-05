#!/usr/bin/env python3

import os, glob

if os.system('git diff --stat --exit-code .'):
    print('There are uncommitted changes in this directory.  Commit these before formatting!')
    exit(0)

for f in glob.iglob('**/*.cpp', recursive=True):
    print('formatting', f)
    os.system('astyle --style=google --indent=spaces=2 %s' % f)
