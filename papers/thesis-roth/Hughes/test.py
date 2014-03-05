from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import Hughes as H
import minmax
import plotHughes as pH

T = 610
prefac = 0.360262397181

pH.plotphi(T,prefac)
plt.show()

# print minmax.maximize(H.Phi,T,prefac,prefac)

# def foo(a, x, c):
#     return x**2

# print minmax.minimize(foo, 0, -40, -39, 0)

