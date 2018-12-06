import numpy as np

import matplotlib.pyplot as plt

data = np.loadtxt('FE_vs_T_at_n_plot_hand_data')
T = data[:,0]

fe = data[:,2]
gw = data[:,3]
plt.plot(T, fe)
plt.axhline(0)
plt.xlabel('T')
plt.ylabel('Free energy difference')

plt.figure()
plt.plot(T, gw)
plt.xlabel('T')
plt.ylabel('width')

plt.show()
