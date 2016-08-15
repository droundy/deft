import sys
import os
import pylab as plt
import time


temp=plt.linspace(.6,1.28,20)
temp=[temp[19]]


ID = int(sys.argv[0].split('ID')[1].split('.py')[0])

print ID
print temp

for i in range(0,len(temp)):
        #if i %4 != ID:
        #        continue
        t=time.time()
        print '%d of %d'%(i,len(temp))
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        os.system('python RG_fn.py %f'%temp[i])
        elapsed = time.time()-t
        print(elapsed)
