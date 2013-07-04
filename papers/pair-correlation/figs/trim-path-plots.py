from pylab import *

# Trims the path data to be used for git storage

for eta in [0.1, 0.2, 0.3, 0.4]:
  data = loadtxt("figs/mc/wallsMC-pair-%02.1f-path.dat" %eta)
  trimmed_data = data[900:1800,0:2]
  savetxt("figs/mc/wallsMC-pair-%02.1f-path-trimmed.dat" %eta, trimmed_data, fmt='%g')

