#!/usr/bin/python2

from __future__ import division
import os, numpy, sys

assert(not os.system("fac ../../../liquid-vapor-monte-carlo"))

i = eval(sys.argv[1])
#RG recursion level

ww = float(sys.argv[2])
#arg ww = [1.3, 1.5, 2.0, 3.0]

L = float((sys.argv[3]))
#arg L = [5.0]

L_i = L*(2**i)
#incorperation of RG parameter i

Ns = eval(sys.argv[4])
#arg Ns = [range(2,10)]

print("\n----------------------LIQUID VAPOR MONTE CARLO------------------------------\n\n")

Overwrite = False
if '-O' in sys.argv:  # Check for overwrite flag in arguments
    Overwrite = True

have_srun = '--srun' in sys.argv
max_non_srun_jobs = 4
jobs_submitted = 0

for N in Ns:
    dirname = 'scrunched-ww%4.2f-L%04.2f/i%01d/N%03d' % (ww,L,i,N)
    if Overwrite:
        os.system('rm -rf ' + dirname)
    print('mkdir -p ' + dirname)
    os.system('mkdir -p '+ dirname)
    #filename = 'ww%4.2f-L%04.2f-N%03d' % ( ww, L, N)
    filename = 'lv-data'
    iterations = 1000000

    output_file_path = dirname+'/'+filename+'-E.dat'
    if not os.path.isfile(output_file_path):
        print "Was checking for", output_file_path
        if have_srun:
            cmd = 'srun -J lv-ww%4.2f-L%04.2f/i%01d/N%03d' % (ww,L,i,N)
        else:
            cmd = ''
        cmd += ' nice -19' # don't hog the CPU
        cmd += ' ../../../liquid-vapor-monte-carlo'
        cmd += ' --filename %s' % filename
        cmd += ' --N %d' % N
        cmd += ' --min-T 0.5'
        cmd += ' --dir %s' % dirname
        cmd += ' --min-samples 10000'
        cmd += ' --tmi'
        cmd += ' --ww %g' % ww
        cmd += ' --iterations %d' %iterations
        cmd += ' --lenx %g --leny %g --lenz %g' % (L_i,L_i,L_i)
        cmd += ' > %s/%s.out 2>&1 &' % (dirname, filename)
        # cmd = ("../../../liquid-vapor-monte-carlo --filename %s --N %d --min_T 0.1 --min_samples 10000 --golden --ww %g --iterations %d --lenx %g --leny %g --lenz %g --data_dir .  > %s.out 2>&1 &" %
        #        (filename, N, ww, iterations, L, L, L, filename))
        print(cmd)
        os.system(cmd)
        jobs_submitted += 1
        if not have_srun and jobs_submitted >= max_non_srun_jobs:
            print("I think this is enough run-monte-carlo processes for now.  No fork bombs!!!\n")
            exit(0)
    else:
        print "You're trying to overwrite: \n %s \n Use the flag -O in order to do so.\n"% output_file_path

if jobs_submitted == 0:
    print "Nothing needs to be done from run-monte-carlo!"
elif jobs_submitted == 1:
    print "I submitted 1 run-monte-carlo job."
else:
    print "I submitted %d run-monte-carlo jobs." % jobs_submitted
