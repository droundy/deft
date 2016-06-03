import sys
import os
import pylab as plt
import time

if len(sys.argv) != 3:
        print("usage: python %s ID num_IDs" % sys.argv[0])
        exit(1)

temp=plt.linspace(0.6,1.28,20)

ID = int(sys.argv[1])
num_IDs = int(sys.argv[2])

print ID
print temp

for i in range(0,len(temp)):
        if i % num_IDs != ID:
                continue
        t=time.time()
        print '%d of %d'%(i,len(temp))
        for j in range(5):
                assert not os.system('python RG_fn.py %f'%temp[i])
        elapsed = time.time()-t
        print(elapsed)
