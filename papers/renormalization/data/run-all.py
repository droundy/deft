#!/usr/bin/python

import sys, os

if len(sys.argv) < 4:
    print "usage: %s RECURSION-LEVEL WELL-WIDTH LENGTH" % (sys.argv[0])
    exit(1)

# Generate a script named rerun.sh, which can be used to run this
# script again with the same arguments.  This is totally gratuitous,
# but can be convenient when I forget which arguments to use.
with open('rerun.sh', 'w') as f:
    f.write('python2 "' + '" "'.join(sys.argv) + '"\n')

i = eval(sys.argv[1])
#RG recursion level

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float((sys.argv[3]))
#arg L = [5.0]

if len(sys.argv) >= 5 and sys.argv[4][0] == '[': # Hokey!!!
    Ns = eval(sys.argv[4])
    #arg Ns = [range(2,10)]
else:
    Ns = range(2,4*(8**i)+1)

have_srun = os.system('squeue') == 0

if '-O' in sys.argv:  # Check for overwrite flag in arguments
    cmd_mc = "python run-monte-carlo.py %d %g %g '%s' -O" %(i,ww,L,Ns)
    if have_srun:
        cmd_mc += ' --srun'
    print(cmd_mc)
    os.system(cmd_mc)

    cmd_fe = 'python run-absolute.py %d %g %g "%s" -O' % (i,ww,L,Ns)
    if have_srun:
        cmd_fe += ' --srun'
    print(cmd_fe)
    os.system(cmd_fe)
else:
    if '--no-abs' not in sys.argv:
        cmd_fe = 'python run-absolute.py %d %g %g "%s"' % (i,ww,L,Ns)
        if have_srun:
            cmd_fe += ' --srun'
        print(cmd_fe)
        os.system(cmd_fe)

    cmd_mc = "python run-monte-carlo.py %d %g %g '%s'" %(i,ww,L,Ns)
    if have_srun:
        cmd_mc += ' --srun'
    print(cmd_mc)
    os.system(cmd_mc)

