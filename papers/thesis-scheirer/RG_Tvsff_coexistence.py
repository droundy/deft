import sys
import os
import pylab as plt
import time





def generate_data():
        temp=plt.linspace(0.6,1.28,20)
        
        for i in range(4,7):
                os.system('python smooth_test_V4.py %0.20f %f'%(temp[i],i))
                print temp[i]
                print i

        

        



generate_data()
