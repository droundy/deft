from __future__ import division
import Hughes as H

# Author: Dan Roth
# E-mail: Daniel.Edward.Roth@gmail.com
# Date: Nov 2013

# Running plotHughes.py will display plots of H.P(T,n) and H.Phi(T,n)
# Looking at the Phi plot we can see what our initial guess should be
# NOTE: the plot is in cgs units and Hughes.py calculates in atomic
# Divide a reading on the plot by H.conv_n to get into atomic units for the calculation

def minimize(func,T,a,c,prefac):

    """Find the minimum of a function in the range [a,b]; return a tuple (abscissa, ordinate).

    func is H.Phi, a function of T, n and prefac
    T is the temp in Kelvin
    a and c are the initial guesses of density (make sure to put into atomic units!)
    prefac is the prefactor term (determines chemical potential)

    """

    fa = func(T,a,prefac)
    fc = func(T,c,prefac)

    b = 0.6*a + 0.4*c
    fb = func(T,b,prefac)
    while fb > fa:
        c = b
        fc = fb
        b = 0.6*a + 0.4*c
        fb = func(T,b,prefac)

    # golden section search
    n0 = a
    n3 = c
    if abs(c - b) > abs(b - a):
        n1 = b
        n2 = b + 0.4*(c - b)
    else:
        n1 = b - 0.4*(b - a)
        n2 = b

    tol = 1.e-6
    f1 = func(T,n1,prefac)
    f2 = func(T,n2,prefac)

    while abs(n3-n0) > tol*(abs(n1)+abs(n2)):
        if f2 < f1: # Shift everything over one way
            n0 = n1
            n1 = n2
            n2 = 0.4*n3 + 0.6*n2

            # re-define f_i to be accurate to new n_i
            f1 = f2
            f2 = func(T,n2,prefac)
        else: # Shift everything over the other way
            n3 = n2
            n2 = n1
            n1 = 0.4*n0 + 0.6*n1

            # re-define f_i to be accurate to the new n_i
            f2 = f1
            f1 = func(T,n1,prefac)

    if f1 < f2:
        nmin = n1
        fmin = f1
    else:
        nmin = n2
        fmin = f2

    return nmin,fmin

def maximize(func,T,a,c,prefac):
    def negfunc(T,n,x):
        return -func(T,n,x)
    n,phi = minimize(negfunc,T,a,c,prefac)
    return n

# def maximize(func,T,g,prefac):
#     """Find maximum of function; return the abscissa.

#     func is H.Phi, a function of T,n and prefac
#     T is the temp in Kelvin
#     g is the initial guess; make sure it's in atomic units
#     prefac is the prefactor term

#     """

#     tol = 1e-8
#     delta = 1e-3
#     Delta = delta*g

#     while func(T,g+Delta,prefac)>func(T,g+tol,prefac) or func(T,g+Delta,prefac)<func(T,g-tol,prefac):
#         if func(T,g+Delta,prefac) > func(T,g,prefac):
#             Delta = delta*g
#             g = g + Delta

#         elif func(T,g+Delta,prefac) < func(T,g,prefac):
#             Delta = delta*g
#             g = g - Delta

#     return g
