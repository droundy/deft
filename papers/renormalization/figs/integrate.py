from __future__ import division
import numpy as np

# Program to do simple integrations

# Author: Dan Roth
# Email: daniel.edward.roth@gmail.com
# Date: January 2014

# very simple
def simple(func,a,b):
    # func should be a function of x

    # initialization
    integral = 0

    # step size
    dx = 1e-3

    # integrate
    for x in np.arange(a,b,dx):
        integral += func(x)*dx

    return integral

# Trapezoidal rule
def trapezoid(func,a,b,n):
    '''func is a function of x
    [a,b] are the limits of integration
    n is the number of segments desired'''

    h = (b-a)/n
    s = func(a) + func(b)

    for k in range (1,n,1):
        s += 2*func(a + k*h)

    return h/2*s


# Midpoint rule
def midpoint(func,a,b,n):
    ''' func is a function of x
    [a,b] are the limits of integration
    n is the number of segments desired '''
    dx = (b-a)/n
    s = 0
    for i in xrange(n):
        s += func(dx*(i + .5))
    return s*dx


# Simpson's rule
def simpson(func,a,b,n):
    '''func is a function of x
    [a,b] are the limits of integration
    n is the number of segments desired
    n must be evnen'''

    h = (b - a)/n
    s = func(a) + func(b)

    for i in range (1,int(n/2),1):
        s += func(a + (2*i-2)*h) + 4*func(a + (2*i-1)*h) + func(a + 2*i*h)

    return h/3*s

def testfunc(x):
    return np.exp(x)

def testfunc_antideriv(x):
    return np.exp(x)

def anal(func,a,b):
    return testfunc(b) - testfunc(a)

# n = 1000
# print 'test e^x from [0,2],%d points'%n
# print '  simple:',simple(testfunc,0,2)
# errsimple = (simple(testfunc,0,2) - anal(testfunc,0,2))/anal(testfunc,0,2)
# print '    err:', errsimple
# print '  trapezoid:',trapezoid(testfunc,0,2,n)
# errtrap = (trapezoid(testfunc,0,2,n) - anal(testfunc,0,2))/anal(testfunc,0,2)
# print '    err:', errtrap
# print '  simpson:',simpson(testfunc,0,2,n)
# errsimpson = (simpson(testfunc,0,2,n) - anal(testfunc,0,2))/anal(testfunc,0,2)
# print '    err:',errsimpson
# print '  anal:',anal(testfunc,0,2)
