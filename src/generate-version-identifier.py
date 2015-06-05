#!/usr/bin/python3

import subprocess, os, sys

fname = sys.argv[1]

try:
    name = subprocess.check_output(['git','describe','--dirty']).decode(encoding='UTF-8')[:-1]
except subprocess.CalledProcessError:
    name = 'unknown'

# We write to a temporary file that we then rename, so that there will
# never be a point in time where fname does not have valid contents.
# This does assume posix file system semantics, and will thus fail on
# Windows.  But so will everything else.

f = open(fname+'~', 'w')
f.write('''static inline const char *version_identifier() {
   return "%s";
}
''' % name)
f.close()

os.rename(fname+'~', fname)
