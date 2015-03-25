#!/usr/bin/python3

import subprocess, os

try:
    name = subprocess.check_output(['git','describe','--dirty']).decode(encoding='UTF-8')[:-1]
except subprocess.CalledProcessError:
    name = 'unknown'

print('''static inline const char *version_identifier() {
   return "%s";
}''' % name)
