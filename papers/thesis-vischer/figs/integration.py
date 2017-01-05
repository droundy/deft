import matplotlib
import numpy as np
matplotlib.use('Agg')
import matplotlib.pyplot as plt


matplotlib.rc('font', **{'family': 'serif', 'serif':['Computer Modern'], 'weight':'bold', 'size': '32'} )
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['legend.fontsize']=28

Ts_ti = [8, 8]
etas_ti = [0, 6]

Ts_mc = [8, 3]
etas_mc = [6,6]

xticks = [0,6]
yticks = [8,3]
T_ticks = ['$T_0=\infty$', '$T_f$']
eta_ticks =['$\eta_0 = 0$', '$\eta_f$']


fig = plt.figure('F')

plt.plot(etas_ti, Ts_ti, "k-o",label='TI', linewidth=4.0)
plt.plot(etas_mc, Ts_mc, "r-o", label = 'SW-MC', linewidth=4.0)

#plt.scatter(etas_ti,Ts_ti)
#plt.scatter(etas_mc,Ts_mc)


#ax = fig.add_subplot(111)
#for axis in ['top', 'bottom','left','right']:
#    ax.spines[axis].set_linewidth(2)


#ax=fig.add_subplot(111)
#ax.text(9,3,"$\eta_0 \to \eta_f$",fontsize= 28)
#ax.text(1.25,1.75,u'',fontsize=28 )

plt.title('Calculating $F_{exc}$ for square well liquid', fontsize = 24, fontweight='bold')
plt.tight_layout(pad=1.5)
plt.legend(loc='best')
plt.xticks(xticks,eta_ticks)
plt.xlim(-.5,10)
plt.ylim(0,10)
plt.yticks(yticks, T_ticks)
plt.ylabel('$T$')
plt.xlabel('$\eta$')
plt.tick_params(axis='both', labelsize=20)

plt.savefig('figs/integration.pdf')
