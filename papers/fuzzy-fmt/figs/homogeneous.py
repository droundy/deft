#!/usr/bin/python

from __future__ import division
# We need the following two lines in order for matplotlib to work
# without access to an X server.
import matplotlib, sys
if 'show' not in sys.argv:
  matplotlib.use('Agg')
from pylab import *
from scipy.special import erf
import os, glob
import styles

matplotlib.rc('text', usetex=True)

#Constants and variables
#k_b = 8.6173324*10**(-5) # in eV
dT = .001
#Temp_max = 600 #in Kelvin
#Temp = arange(.001, .1 + dT/2, dT)
epsilon = 1
R = 1# in Angstroms
density = arange(0, .8 - .001/2, .0001)/(4*pi/3)

eta = density*4*pi/3
#P_cs = density*.001*(1+eta+eta**2-eta**3)/(1-eta)**3
P_cs = density*(1+eta+eta**2)/(1-eta)**3
#plot(eta,P_cs/.001, 'k',linewidth=2, label = 'Hard spheres')

Temp = 0.00001
eta = density*4*pi/3
#P_cs = density*(1+eta+eta**2-eta**3)/(1-eta)**3
P_cs = density*(1+eta+eta**2)/(1-eta)**3
P_cs = density

bhdata = loadtxt('figs/homogeneous-bh.dat')
bhtemp = bhdata[0,1:]
bhnred = bhdata[1:,0]
bhpressures = bhdata[1:, 1:]
bhnans = bhpressures != bhpressures
bhpressures[bhnans] = 1e30

erfdata = loadtxt('figs/homogeneous.dat')
erftemp = erfdata[0,1:]
erfnred = erfdata[1:,0]

pressures = erfdata[1:, 1:]
mynans = pressures != pressures
pressures[mynans] = 1e30

temperatures, n_reduced = meshgrid(erftemp, erfnred)

bhtemperatures, bhn_reduced = meshgrid(bhtemp, bhnred)

# for j in arange(1,len(erfdata[0,:])):
#   erfpressure = erfdata[:,j]
#   plot(erfnred[1:], erfpressure[1:], styles.color[erftemp[j]]+'-', label='$kT/\epsilon$=%g' % erftemp[j])

dp = 0.25
levels = arange(dp, 2 + dp/2, dp)

CS = contour(temperatures, n_reduced, pressures, levels)
clabel(CS, inline=1, fontsize=10)

title('FIXME: Check on weighting functions for homogeneous at low $T$')

CS = contour(temperatures, n_reduced, n_reduced*temperatures, levels, linestyles='dashed')
clabel(CS, inline=1, fontsize=10)

ylabel('n')
xlabel('T')
ylim(ymax=1.0)

savefig('figs/n-vs-T.pdf', bbox_inches=0)

def label_line(line, label_text, near_i=None, near_x=None, near_y=None, rotation_offset=0, offset=(0,0)):
    """call 
        l, = plt.loglog(x, y)
        label_line(l, "text", near_x=0.32)
    """
    def put_label(i):
        """put label at given index"""
        i = min(i, len(x)-2)
        dx = sx[i+1] - sx[i]
        dy = sy[i+1] - sy[i]
        rotation = np.rad2deg(math.atan2(dy, dx)) + rotation_offset
        pos = [(x[i] + x[i+1])/2. + offset[0], (y[i] + y[i+1])/2 + offset[1]]
        plt.text(pos[0], pos[1], label_text, size=9, rotation=rotation, color = line.get_color(),
        ha="center", va="center", bbox = dict(ec='1',fc='1'))

    x = line.get_xdata()
    y = line.get_ydata()
    ax = line.get_axes()
    if ax.get_xscale() == 'log':
        sx = np.log10(x)    # screen space
    else:
        sx = x
    if ax.get_yscale() == 'log':
        sy = np.log10(y)
    else:
        sy = y/50.0

    # find index
    if near_i is not None:
        i = near_i
        if i < 0: # sanitize negative i
            i = len(x) + i
        put_label(i)
    elif near_x is not None:
        for i in range(len(x)-2):
            if (x[i] < near_x and x[i+1] >= near_x) or (x[i+1] < near_x and x[i] >= near_x):
                put_label(i)
    elif near_y is not None:
        for i in range(len(y)-2):
            if (y[i] < near_y and y[i+1] >= near_y) or (y[i+1] < near_y and y[i] >= near_y):
                put_label(i)
    else:
        raise ValueError("Need one of near_i, near_x, near_y")

figure()
skipby = 50 # number of temperatures to skip
xlim(xmin=0.,xmax=1.2)
ylim(ymin=0.001, ymax = 30)
near_x = 0.85
for i in range(skipby-1,len(n_reduced[0,:]),skipby):
  line, = semilogy(n_reduced[:,i], pressures[:,i], label='T=%g' % temperatures[0,i])
  temp_name = '$T = %g$' % temperatures[0,i]
  if temperatures[0,i] < 1 or not ('.' in temp_name):
    label_line(line, temp_name, near_x=near_x)
    if temperatures[0,i] < 1:
      near_x += 0.05
    else:
      near_x += 0.01
  #plot(1/n_reduced[:,i], temperatures[0,i]*n_reduced[:,i], ':')

xlabel('$1/n$')
ylabel('$p$')

savefig('figs/p-vs-V.pdf', bbox_inches=0)

figure()
skipby = 10 # number of temperatures to skip
xlim(xmin=0.,xmax=1.2)
ylim(ymin=0, ymax = 15)
near_x = 0.85
for i in range(skipby-1,len(n_reduced[0,:]),skipby):
  line, = plot(n_reduced[:,i], pressures[:,i], label='T=%g' % temperatures[0,i])
  temp_name = '$T = %g$' % temperatures[0,i]
  if temperatures[0,i] < 1 or not ('.' in temp_name):
    label_line(line, temp_name, near_x=near_x)
    if temperatures[0,i] < 1:
      near_x += 0.05
    else:
      near_x += 0.01
  #plot(1/n_reduced[:,i], temperatures[0,i]*n_reduced[:,i], ':')

xlabel('$1/n$')
ylabel('$p$')

savefig('figs/p-vs-n.pdf', bbox_inches=0)

figure()

Nlabel = []
Tlabel = []
Plabel = []
MCPlabel = []

for i in range(9,len(n_reduced[:,0]),10)[:10]:
  plot(temperatures[i,:], pressures[i,:],
       styles.density_color(n_reduced[i,0])+'-')
  plot(bhtemperatures[i,:], bhpressures[i,:],
       styles.density_color(n_reduced[i,0])+':')
  if i > 9:
    Nlabel.append(n_reduced[i,0])
    Tlabel.append(2.75 - 2.7*(n_reduced[i,0] - 0.2)**0.9)
    Plabel.append(pressures[i,:][temperatures[i,:] >= Tlabel[-1]][0])
    MCPlabel.append(0.0)

Nlabel = np.array(Nlabel)
Tlabel = np.array(Tlabel)
Plabel = np.array(Plabel)
MCPlabel = np.array(MCPlabel)

# input: 'papers/fuzzy-fmt/figs/mcwca-*-*.dat.prs'
# for rd in arange(0.70, 0.35, -0.1):
#   temps = []
#   pressures = []
#   print 'could use:', glob.glob('figs/mcwca-%.4f-*.dat.prs' % (rd))
#   for temp in arange(0.1, 3.05, 0.1):
#     fname = 'figs/mcwca-%.4f-%.4f.dat.prs' % (rd, temp)
#     if os.path.exists(fname):
#       temps.append(temp)
#       p = loadtxt(fname)
#       # we account below for the fact that the MC code does not spit
#       # out a reduced pressure, so we need to do a conversion.
#       pressures.append(p*2.0**(5/2.0))
#   if len(temps) > 0:
#     plot(temps, pressures, styles.density_color(rd)+'--',
#          label='MC $n^{*} = %g$' % rd)

# input: 'papers/fuzzy-fmt/figs/mcfcc-1.0000-1.5000.dat.prs'
# input: 'papers/fuzzy-fmt/figs/mcfcc-1.0000-2.0000.dat.prs'
# input: 'papers/fuzzy-fmt/figs/mcfcc-0.9000-2.0000.dat.prs'
# input: 'papers/fuzzy-fmt/figs/mcfcc-0.8000-2.0000.dat.prs'
# input: 'papers/fuzzy-fmt/figs/mcfcc-0.8000-1.0000.dat.prs'
# input: 'papers/fuzzy-fmt/figs/mcfcc-*-*.dat.prs'
for rd in arange(1.00, 0.25, -0.1):
  temps = []
  pressures = []
  # print 'could use:', glob.glob('figs/mcfcc-%.4f-*.dat.prs' % (rd))
  for temp in arange(0.05, 3.05, 0.05):
    fname = 'figs/mcfcc-%.4f-%.4f.dat.prs' % (rd, temp)
    if os.path.exists(fname):
      temps.append(temp)
      p = loadtxt(fname)
      # we account below for the fact that the MC code does not spit
      # out a reduced pressure, so we need to do a conversion.
      pressures.append(p*2.0**(5/2.0))
  if len(temps) > 0:
    plot(temps, pressures, styles.density_color(rd)+'--')
    if len(np.argwhere(np.abs(Nlabel - rd) < 0.001)) == 1:
      print 'We have MC for rd', rd
      pressures = np.array(pressures)
      temps = np.array(temps)
      i = np.argwhere(np.abs(Nlabel - rd) < 0.001)[0][0]
      j = np.argwhere(temps > Tlabel[i])[0][0]
      MCPlabel[i] = pressures[j]


for i in range(0,len(Tlabel)):
  if MCPlabel[i] == 0:
    MCPlabel[i] = Plabel[i]
  text(Tlabel[i], 0.5*(Plabel[i] + MCPlabel[i]), '$n^{*} = %g$' % Nlabel[i],
       color=styles.density_color(Nlabel[i]),
       backgroundcolor='w',
       va="center", ha="center")

xlim(xmin=0, xmax = 3)
ylim(ymin=0, ymax = 25)
xlabel('$T^*$')
ylabel('$p^*$')

plot([],[], 'k-', label='SFMT theory')
plot([],[], 'k--', label='Monte Carlo')
plot([],[], 'k:', label='BH')
legend(loc='best', frameon=False)

savefig('figs/p-vs-T.pdf', bbox_inches=0)

show()
