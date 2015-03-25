#!/usr/bin/python3

import subprocess, os

try:
    #name = subprocess.check_output(['git','describe','--dirty']).decode(encoding='UTF-8')[:-1]
    # the above returns an error on quipu
    # while adding try/except at least lets scons build with the above, it won't get us git version info (returning name = 'unknown')
    # we might want to use the following instead, which properly gets git version info on at least two machines (namely quipu and MPerlin's computer)
    # add '--short' before 'HEAD' for a shorter snippet, if desired
    name = subprocess.check_output(['git','rev-parse','HEAD']).decode(encoding='UTF-8')[:-1]
except subprocess.CalledProcessError:
    name = 'unknown'

print('''static inline const char *version_identifier() {
   return "%s";
}''' % name)
