import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt


matplotlib.rc('font', **{'family': 'serif', 'serif':['Computer Modern'], 'weight':'bold', 'size': '32'} )
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['legend.fontsize']=28

rs = [.5,.5,1.5,1.5,5]
rs2 = [.5, 1.5]
vs2 = [1.5, 1.5]
rticks = ['$0$','$\sigma$', '$\lambda \sigma$']
vticks = ['$- \epsilon$', '$0$']

vs = [10,.5,.5,1.5,1.5]
val = 0

fig = plt.figure('s-w')
#plt.scatter(rs,vs)
plt.plot(rs,vs, color="black",label='$\Phi_{SW}$', linewidth=2.0)
plt.plot(rs2, vs2, 'k--', label = '$\Phi_{HS}$', linewidth=3.0)

#ax = fig.add_subplot(111)
#for axis in ['top', 'bottom','left','right']:
#    ax.spines[axis].set_linewidth(2)
    
plt.title('Potential for hard sphere and square well fluids', fontsize = 24, fontweight='bold')
plt.tight_layout(pad=1.5)
plt.legend(loc='best')
plt.xticks([0,.5,1.5],rticks)
plt.xlim(0,3)
plt.ylim(0,4)
plt.yticks([.5,1.5], vticks)
plt.ylabel('$\mathbf{\Phi_{12}}$')
plt.xlabel('$|\mathbf{r_1} - \mathbf{r_2}|$')
plt.tick_params(axis='both', labelsize=20)

plt.savefig('figs/phi12.pdf')
