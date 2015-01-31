#!/usr/bin/python3

print("""#!/usr/bin/python3
""")

import string, os

os.mkdir('testing-flags');
with open('testing-flags/test.c', 'w') as f:
    f.write("""int main() {
  return 0;
}
""")
flags = ''
for flag in ['-Wall', '-Werror', '-O3', '-thisdosntexist', '-std=c++11', '-Isrc', '-Iinclude', '-g']:
    if not os.system('cd testing-flags && g++ %s %s -c test.c' %
                     (flags, flag)):
        flags += ' ' + flag
    else:
        print('# g++ cannot use flag:', flag)
if len(flags) > 0:
    flags = flags[1:]
linkflags = ''
for flag in ['-lpopt', '-lprofiler', '-lfftw3']:
    if not os.system('cd testing-flags && g++ %s %s -o test test.c' %
                     (flags, flag)):
        linkflags += ' ' + flag
    else:
        print('# g++ linking cannot use flag:', flag)
if len(linkflags) > 0:
    linkflags = linkflags[1:]
os.system('rm -rf testing-flags')

print('cxx=', repr('g++')) # in future, we should set this dynamically
print('cxxflags=', repr(flags))
print('linkflags=', repr(linkflags))
