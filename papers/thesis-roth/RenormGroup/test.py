from __future__ import division
import numpy as np
import Hughes
from scipy import integrate

def fbarD(T,nl,x):
    if x == 1:
        return x*10
    else:
        return x/10

# def ID(T,nl,x,i):
#     def func(x):
#         return fbarD(T,nl,x) + UbarD(T,nl,x)
#     def integrand(x):
#         return np.exp(func(x))
#     return integrate.quad(integrand,0,nl)[0]

def ID(T,nl,x,i):
    def func(x):
        return np.exp(fbarD(T,nl,x) + UbarD(T,nl,x))
    return integrate.quad(func,0,nl)[0]

def UbarD(T,nl,x):
    return x

def func(x):
    result = 0
    for i in range(x+1):
        result += i
        print result

func(5)
