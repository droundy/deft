import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt


matplotlib.rc('font', **{'family': 'serif', 'serif':['Computer Modern'], 'weight':'bold', 'size': '32'} )
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['xtick.labelsize'] = 50
matplotlib.rcParams['ytick.labelsize'] = 50

c= [1.5,1.5]
t= [1,1]

l1x = np.arange(0,t[0],.01)
l1y = map(lambda x: x**2, l1x)

l2x = np.arange(t[0],t[0]+.5,.01)
l2y = map(lambda x: x*20 -20+t[0], l2x)

l3x = np.arange(t[0],c[0],.01)
l3y = map(lambda x: 2*(x-t[0])**2 +t[1], l3x)


xticks = ['$0$','$T_t$', '$T_c$']
yticks = ['$P_t$', '$P_c$']



fig=plt.figure('phase')
plt.scatter([t[0],c[0]],[t[1],c[1]])
plt.plot(l1x,l1y, color="black", linewidth=4)
plt.plot(l2x,l2y, color="black", linewidth=4)
plt.plot(l3x,l3y, color="black", linewidth=4)
plt.plot([t[0],t[0]],[t[1],0], 'k--', lw=2)
plt.plot([t[0],0],[t[1],t[1]], 'k--', lw=2)
plt.plot([c[0],c[0]],[c[1],0], 'k--', lw=2)
plt.plot([c[0],0],[c[1],c[1]],'k--',lw=2)

ax=fig.add_subplot(111)
ax.text(.5,2.5,u'S',fontsize= 28)
ax.text(1.25,1.75,u'L',fontsize=28 )
ax.text(1.75,.75,u'G',fontsize=28 )
ax.text(2,2,u'F', fontsize=28)



plt.tight_layout(pad=1)
plt.xticks([0,t[0],c[0]],xticks, fontsize=50)
plt.xlim(0,3)
plt.ylim(0,4)
plt.yticks([t[1],c[1]], yticks, fontsize=50)
plt.ylabel('Pressure')
plt.xlabel('Temperature')
plt.tick_params(axis='both', labelsize=20)

plt.savefig('figs/phase-new.pdf')
