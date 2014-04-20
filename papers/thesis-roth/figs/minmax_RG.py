from __future__ import division

import numpy as np

# Author: Dan Roth
# E-mail: Daniel.Edward.Roth@gmail.com
# Date: Nov 2013

def minimize(func,T,a,c,npart,i):

    """Find the minimum of a function in the range [a,c]; return a tuple (abscissa, ordinate).

    func is grand free energy, a function of T, n and some nparticular (for calculating mu)
    T is the temp in Kelvin
    a and c are the initial guesses of density (make sure to put into atomic units!)
    npart is the term that determines chemical potential
    i is the iteration depth

    """

    fa = func(T,a,npart,i)
    fc = func(T,c,npart,i)

    b = 0.6*a + 0.4*c
    b_old = a
    fb = func(T,b,npart,i)
    dx = (c - a)/100000

    if a == c:
        return a, fa
    possible_bs = np.arange(a+dx, c-dx, dx)
    f_possible_bs = func(T, possible_bs, npart,i)
    for bindex in xrange(1, len(possible_bs)):
        if f_possible_bs[bindex] < fb:
            b = possible_bs[bindex]
            fb = f_possible_bs[bindex]
    if fb > fa:
        return a, fa
    if fb > fc:
        return c, fc

    # golden section search
    n0 = a
    n3 = c
    if abs(c - b) > abs(b - a):
        n1 = b
        n2 = b + 0.4*(c - b)
    else:
        n1 = b - 0.4*(b - a)
        n2 = b

    tol = 1.e-10
    f1 = func(T,n1,npart,i)
    f2 = func(T,n2,npart,i)

    while abs(n3-n0) > tol*(abs(n1)+abs(n2)):
        if f2 < f1: # Shift everything over one way
            n0 = n1
            n1 = n2
            n2 = 0.4*n3 + 0.6*n2

            # re-define f_i to be accurate to new n_i
            f1 = f2
            f2 = func(T,n2,npart,i)
        else: # Shift everything over the other way
            n3 = n2
            n2 = n1
            n1 = 0.4*n0 + 0.6*n1

            # re-define f_i to be accurate to the new n_i
            f2 = f1
            f1 = func(T,n1,npart,i)

    if f1 < f2:
        nmin = n1
        fmin = f1
    else:
        nmin = n2
        fmin = f2

    return nmin,fmin

def maximize(func,T,a,c,npart,i):
    def negfunc(T,n,x,ineg):
        return -func(T,n,x,ineg)
    n,phi = minimize(negfunc,T,a,c,npart,i)
    return n
