#!/usr/bin/python3

import subprocess, os

if os.path.exists('.git/HEAD'):
    name = subprocess.check_output(['git', 'describe', '--dirty']).decode(encoding='UTF-8')[:-1]
else:
    name = 'unknown'

print('''static inline const char *version_identifier() {
   return "%s";
}''' % name)
