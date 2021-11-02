#!/bin/sh

set -ev

fac new-radial-lj.mkdat new-radial-wca.mkdat # new-bh-radial-wca.mkdat new-bh-radial-lj.mkdat

# ../../../src/new-radial-lj.mkdat 

# ./new-radial-wca.mkdat DENSITY TEMPERATURE 
rq run -J radial-6-10 --mem 88 -o new-data/radial-wca-0.6000-10.00.out ./new-radial-wca.mkdat 0.6 10.0 
rq run -J radial-6-5 --mem 88 -o new-data/radial-wca-0.6000-5.00.out ./new-radial-wca.mkdat 0.6 5.0
rq run -J radial-6-2.5 --mem 88 -o new-data/radial-wca-0.6000-2.50.out ./new-radial-wca.mkdat 0.6 2.5
rq run -J radial-6-1 --mem 88 -o new-data/radial-wca-0.6000-1.00.out ./new-radial-wca.mkdat 0.6 1
rq run -J radial-6-0.1 --mem 88 -o new-data/radial-wca-0.6000-0.10.out ./new-radial-wca.mkdat 0.6 0.1

rq run -J radial-1-10 --mem 88 -o new-data/radial-wca-1.0000-10.00.out ./new-radial-wca.mkdat 1 10.0
rq run -J radial-1-5 --mem 88 -o new-data/radial-wca-1.0000-5.00.out ./new-radial-wca.mkdat 1 5.0
rq run -J radial-1-2.5 --mem 88 -o new-data/radial-wca-1.0000-2.50.out ./new-radial-wca.mkdat 1 2.5

#Run new-bh-radial-wca.cpp 		
rq run -J radialbh-6-10 --mem 88 -o new-data/radial-bh-wca-0.6000-10.00.out ./new-bh-radial-wca.mkdat 0.6 10.0 
rq run -J radialbh-6-5 --mem 88 -o new-data/radial-bh-wca-0.6000-5.00.out ./new-bh-radial-wca.mkdat 0.6 5.0     
rq run -J radialbh-6-2.5 --mem 88 -o new-data/radial-bh-wca-0.6000-2.50.out ./new-bh-radial-wca.mkdat 0.6 2.5  
rq run -J radialbh-6-1 --mem 88 -o new-data/radial-bh-wca-0.6000-1.00.out ./new-bh-radial-wca.mkdat 0.6 1     
#(rq run -J radialbh-6-0.1 --mem 88 -o new-data/radial-bh-wca-0.6000-0.10.out ./new-bh-radial-wca.mkdat 0.6 0.1 	NO?)

rq run -J radialbh-1-10 --mem 88 -o new-data/radial-bh-wca-1.0000-10.00.out ./new-bh-radial-wca.mkdat 1 10.0	
rq run -J radialbh-1-5 --mem 88 -o new-data/radial-bh-wca-1.0000-5.00.out ./new-bh-radial-wca.mkdat 1 5.0	
rq run -J radialbh-1-2.5 --mem 88 -o new-data/radial-bh-wca-1.0000-2.50.out ./new-bh-radial-wca.mkdat 1 2.5  



#Run new-radial-lj.cpp 
rq run -J radiallj-0.83890-0.71 --mem 88 -o new-data/radial-lj-0.71-0.83890.out ./new-radial-lj.mkdat 0.83890 0.71 
#rq run -J radiallj-6-5 --mem 88 -o new-data/radial-lj-2.48-0.957.out ./new-radial-lj.mkdat 0.957 2.48     
rq run -J radiallj-0.5844-1.235 --mem 88 -o new-data/radial-lj-1.235-0.5844.out ./new-radial-lj.mkdat 0.5844 1.235  
rq run -J radiallj-1.095-2.48 --mem 88 -o new-data/radial-lj-2.48-1.095.out ./new-radial-lj.mkdat 1.095 2.48     


#Run new-bh-radial-lj.cpp
rq run -J radialbh-lj-0.83890-0.71 --mem 88 -o new-data/radial-bh-lj-0.71-0.83890.out ./new-bh-radial-lj.mkdat 0.83890 0.71
rq run -J radialbh-lj-1.095-2.48 --mem 88 -o new-data/radial-bh-lj-2.48-1.095.out ./new-bh-radial-lj.mkdat  1.095 2.48     
rq run -J radialbh-lj-0.5844-1.235 --mem 88 -o new-data/radial-bh-lj-1.235-0.5844.out ./new-bh-radial-lj.mkdat 0.5844 1.235  
    



# def runme(reduced_density, temperature):
#     outfilename = figsdir+'/new-data/radial-lj-%06.4f-%04.2f.out' % (temperature, reduced_density)
#     #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
#     os.system("%s %s/new-radial-lj.mkdat %g %g > %s 2>&1 &" %
#            (srun(reduced_density, temperature), figsdir, reduced_density, temperature, outfilename))

# def run_wca(reduced_density, temperature):
#     outfilename = figsdir+'/new-data/radial-wca-%06.4f-%04.2f.out' % (temperature, reduced_density)
#     #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
#     cmd = ("%s %s/new-radial-wca.mkdat %g %g > %s 2>&1 &" %
#            (srun(reduced_density, temperature).replace('name', 'new-wca'), figsdir, reduced_density, temperature, outfilename))
#     print(cmd)
#     os.system(cmd)

# def run_bh_wca(reduced_density, temperature):
#     outfilename = figsdir+'/new-data/radial-bh-wca-%06.4f-%04.2f.out' % (temperature, reduced_density)
#     #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
#     cmd = ("%s %s/new-bh-radial-wca.mkdat %g %g > %s 2>&1 &" %
#            (srun(reduced_density, temperature).replace('name', 'new-bh-wca'), figsdir, reduced_density, temperature, outfilename))
#     print(cmd)
#     os.system(cmd)

# def run_bh_lj(reduced_density, temperature):
#     outfilename = figsdir+'/new-data/radial-bh-lj-%06.4f-%04.2f.out' % (temperature, reduced_density)
#     #system("%s %s/soft-sphere.mkdat %g %g > %s 2>&1 &" %
#     cmd = ("%s %s/new-bh-radial-lj.mkdat %g %g > %s 2>&1 &" %
#            (srun(reduced_density, temperature).replace('name', 'new-bh-lj'), figsdir, reduced_density, temperature, outfilename))
#     print(cmd)
#     os.system(cmd)


# runme(0.83890, 0.71)

# runme(0.957, 2.48)

# runme(0.5844, 1.235)

# runme(1.095, 2.48)
# run_wca(0.6, 10.0)
# run_wca(0.6, 5.0)
# run_wca(0.6, 2.5)

#for temp in [2.5]:
#    run_bh_wca(0.6, temp)
#    run_bh_wca(1.0, temp)
    
#run_bh_wca(0.6, 2.5)
#run_bh_wca(0.6, 0.1)

# run_bh_lj(0.83890, 0.71)
# run_bh_lj(1.095, 2.48)
# run_bh_lj(0.5844, 1.235)
