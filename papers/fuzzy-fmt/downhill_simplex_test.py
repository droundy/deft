#!/usr/bin/python2
#run this from fuzzy-fmt

import os
import numpy as np
import matplotlib.pyplot as plt
import sys

testdata = np.loadtxt("test_simplex_data_file")
print testdata

plt.title("simplex")
plt.xlabel("fv")
plt.ylabel("gw")
plt.scatter(testdata[0::1,0::2], testdata[0::1,1::2], color='blue')
#print testdata[0::1,0::2]
#print testdata[0::1,1::2]
plt.savefig("saveplot")
#plt.show

#print testdata[0,0::2] #delete
#print testdata[0,1::2] #delete

failflag=0

for i in range(0,len(testdata)):
    print
    print
    print "LOOP", i, "-------------"
    print "point 1: ", testdata[i,0],testdata[i,1]
    print "point 2: ", testdata[i,2],testdata[i,3]
    print "point 3: ", testdata[i,4],testdata[i,5]
    plt.scatter(testdata[i,0::2], testdata[i,1::2], color='blue')
    #compute reflected point
    ##reflectedfv=testdata[0,0]+testdata[1,0]-testdata[2,0] #corresponds to new-melting
    ##reflectedgw=testdata[0,1]+testdata[1,1]-testdata[2,1]
    reflectedfv=testdata[i,0]+testdata[i,2]-testdata[i,4]
    reflectedgw=testdata[i,1]+testdata[i,3]-testdata[i,5]
    #plt.scatter(testdata[i,0::2], testdata[i,1::2], color='blue')
    #plt.scatter(reflectedfv, reflectedgw, color='red')  #how add point to plot?
    #plt.scatter(reflectedfv, reflectedgw, c='red')
    print
    print "reflected fv=", reflectedfv
    print "reflected gw=", reflectedgw
    diffgw_pt2_ref=testdata[i,3]-reflectedgw
    print "diffgw_pt2_ref=", diffgw_pt2_ref
    difffv_pt2_ref=testdata[i,2]-reflectedfv
    print "difffv_pt2_ref=", difffv_pt2_ref
    diffgw_pt3_pt1=testdata[i,5]-testdata[i,1]
    print "diffgw_pt3_pt1=", diffgw_pt3_pt1
    difffv_pt3_pt1=testdata[i,4]-testdata[i,0] 
    print "difffv_pt3_pt1=", difffv_pt3_pt1
    if (diffgw_pt2_ref == diffgw_pt3_pt1 and difffv_pt2_ref==difffv_pt3_pt1) :
        print "*reflected PASSED"
    else: 
        print "!reflected FAILED"
        failflag=1
    
    #compute extended point
    ##extendedfv=((3/2.0)*(testdata[0,0]+testdata[1,0]))-(2.0*testdata[2,0])
    ##extendedgw=((3/2.0)*(testdata[0,1]+testdata[1,1]))-(2.0*testdata[2,1])
    extendedfv=((3/2.0)*(testdata[i,0]+testdata[i,2]))-(2.0*testdata[i,4])
    extendedgw=((3/2.0)*(testdata[i,1]+testdata[i,3]))-(2.0*testdata[i,5])
    print
    print "extended fv=", extendedfv
    print "extended gw=", extendedgw
    diffgw_ext_ref=extendedgw - reflectedgw
    print "diffgw_ext_ref=", diffgw_ext_ref
    difffv_ext_ref=extendedfv - reflectedfv
    print "difffv_ext_ref=", difffv_ext_ref
    diffgw_ref_pt1=reflectedgw-testdata[i,5]
    print "diffgw_ref_pt1=", diffgw_ref_pt1    
    difffv_ref_pt1=reflectedfv-testdata[i,4] 
    print "difffv_ref_pt1=", difffv_ref_pt1
    if (diffgw_ext_ref==(1.0/2)*diffgw_ref_pt1 and difffv_ext_ref==(1.0/2)*difffv_ref_pt1) :
        print "*extended PASSED"
    else:
        print "!extended FAILED"
        failflag=2

    #compute contacted_out point
    ##contracted_out_fv=((3/4.0)*(testdata[0,0]+testdata[1,0]))-((1/2.0)*testdata[2,0])
    ##contracted_out_gw=((3/4.0)*(testdata[0,1]+testdata[1,1]))-((1/2.0)*testdata[2,1])
    contracted_out_fv=((3/4.0)*(testdata[i,0]+testdata[i,2]))-((1/2.0)*testdata[i,4])
    contracted_out_gw=((3/4.0)*(testdata[i,1]+testdata[i,3]))-((1/2.0)*testdata[i,5])
    print
    print "contracted_out fv=", contracted_out_fv
    print "contracted_out gw=", contracted_out_gw
    diffgw_ref_conout=reflectedgw-contracted_out_gw
    print "diffgw_ref_conout=", diffgw_ref_conout
    difffv_ref_conout=reflectedfv-contracted_out_fv
    print "difffv_ref_conout=", difffv_ref_conout
    diffgw_conout_pt3=contracted_out_gw-testdata[i,5]
    print "diffgw_conout_pt3=", diffgw_conout_pt3
    difffv_conout_pt3=contracted_out_fv-testdata[i,4] 
    print "difffv_conout_pt3=", difffv_conout_pt3
    if (diffgw_ref_conout == (1.0/3)*diffgw_conout_pt3 and difffv_ref_conout==(1.0/3)*difffv_conout_pt3) :
        print "*contracted_out PASSED"
    else:
        print "!contracted_out FAILED"
        failflag=3

    #compute contacted_in point
    ##contracted_in_fv=((1/4.0)*(testdata[0,0]+testdata[1,0]))-((1/2.0)*testdata[2,0])
    ##contracted_in_gw=((1/4.0)*(testdata[0,1]+testdata[1,1]))-((1/2.0)*testdata[2,1])
    contracted_in_fv=((1/4.0)*(testdata[i,0]+testdata[i,2]))+((1/2.0)*testdata[i,4])
    contracted_in_gw=((1/4.0)*(testdata[i,1]+testdata[i,3]))+((1/2.0)*testdata[i,5])
    print
    print "contracted_in fv=", contracted_in_fv
    print "contracted_in gw=", contracted_in_gw
    diffgw_ref_conin=reflectedgw-contracted_in_gw
    print "diffgw_ref_conin=", diffgw_ref_conin
    difffv_ref_conin=reflectedfv-contracted_in_fv
    print "difffv_ref_conin=", difffv_ref_conin
    diffgw_conin_pt3=contracted_in_gw-testdata[i,5]
    print "diffgw_conin_pt3=", diffgw_conin_pt3
    difffv_conin_pt3=contracted_in_fv-testdata[i,4] 
    print "difffv_conin_pt3=", difffv_conin_pt3
    if (diffgw_ref_conin == 3*diffgw_conin_pt3 and difffv_ref_conin==3*difffv_conin_pt3) :
        print "*contracted_in PASSED"
    else:
        print "!contracted_in FAILED"
        failflag=4

    #compute shrink_in point
    ##shrink_in_fv=(1/2.0)*(testdata[0,0]+testdata[2,0])
    ##shrink_in_gw=(1/2.0)*(testdata[0,1]+testdata[2,1])
    shrink_in_fv=(1/2.0)*(testdata[i,0]+testdata[i,4])
    shrink_in_gw=(1/2.0)*(testdata[i,1]+testdata[i,5])
    print
    print "shrink_in fv=", shrink_in_fv
    print "shrink_in gw=", shrink_in_gw
    diffgw_pt2_shkin=testdata[i,5]-shrink_in_gw
    print "diffgw_pt2_shkin=", diffgw_pt2_shkin
    difffv_pt2_shkin=testdata[i,4]-shrink_in_fv
    print "difffv_pt2_shkin=", difffv_pt2_shkin
    diffgw_shkin_pt1=shrink_in_gw-testdata[i,1]
    print "diffgw_shkin_pt1=", diffgw_shkin_pt1
    difffv_shkin_pt1=shrink_in_fv-testdata[i,0] 
    print "difffv_shkin_pt1=", difffv_shkin_pt1
    if (diffgw_pt2_shkin == diffgw_shkin_pt1 and difffv_pt2_shkin==difffv_shkin_pt1) :
        print "*shrink_in PASSED"
    else:
        print "!shrink_in FAILED"
        failflag=5

    #compute shrink_out point
    ##shrink_out_fv=(1/2.0)*(testdata[0,0]+testdata[1,0])
    ##shrink_out_gw=(1/2.0)*(testdata[0,1]+testdata[1,1])
    shrink_out_fv=(1/2.0)*(testdata[i,0]+testdata[i,2])
    shrink_out_gw=(1/2.0)*(testdata[i,1]+testdata[i,3])
    print
    print "shrink_out fv=", shrink_out_fv
    print "shrink_out gw=", shrink_out_gw
    diffgw_pt2_shkout=testdata[i,3]-shrink_out_gw
    print "diffgw_pt2_shkout=", diffgw_pt2_shkout
    difffv_pt2_shkout=testdata[i,2]-shrink_out_fv
    print "difffv_pt2_shkout=", difffv_pt2_shkout
    diffgw_shkout_pt1=shrink_out_gw-testdata[i,1]
    print "diffgw_shkout_pt1=", diffgw_shkout_pt1
    difffv_shkout_pt1=shrink_out_fv-testdata[i,0] 
    print "difffv_shkout_pt1=", difffv_shkout_pt1
    if (diffgw_pt2_shkout == diffgw_shkout_pt1 and difffv_pt2_shkout==difffv_shkout_pt1) :
        print "*shrink_out PASSED"
    else: 
        print "!shrink_out FAILED"
        failflag=6
    
    print "---------------------"    
    
  
if (failflag==0):
    print "ALL PASSED"
else:
    print "FAILED"
