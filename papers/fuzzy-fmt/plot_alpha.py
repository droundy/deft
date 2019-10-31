#!/usr/bin/python2

#Run this program from the directory it is listed in
#with command ./plot_alpha.py

from scipy import special
import numpy as np
import matplotlib.pyplot as plt
import math

epsilon=1
sigma=1
PI=np.pi

KbT=np.linspace(0.05, 2000, 2000)

alpha=sigma*np.cbrt(np.sqrt((2/(1+np.sqrt((KbT*np.log(2))/epsilon)))))

#Plot alpha
plt.plot(KbT, alpha, label = 'alpha')
plt.xlabel('KbT')
plt.ylabel('alpha')
plt.title('alpha vs Temp')
plt.legend()

plt.show()   
