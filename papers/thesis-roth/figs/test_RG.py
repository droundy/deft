from __future__ import division
import SW
import RG

T = 5
n = 0.02

print '    SW; T = %d; n = %0.2f; f = %0.4f'%(T,n,SW.f(T,n))
print '\n'

for i in range (2):
    print '    RG; i = %d; T = %d; n = %0.2f; f = %0.4f'%(i,T,n,RG.ftot(T,n,i))
