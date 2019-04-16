#!/usr/bin/python3

#Use this program to find the minimum of a set of data points 
#Run this program from deft/papers/fuzzy-fmt with command python3 test-minimization.py

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt

xmin = 3.0
eps = 2.0
emin = 1.0
sigma = 0.1

def func_to_minimize(x):
    return emin + 4*eps*((xmin/x)**12 - (xmin/x)**6) + eps + np.random.normal()*sigma


def parabola_fit(xs,ys):
    A = np.vstack([xs**2, xs, np.ones(len(xs))]).T
    a, b, c = np.linalg.lstsq(A, ys, rcond=None)[0]
    x0 = -b/(2*a)
    y0 = a*x0**2 + b*x0 + c
    error = np.mean(np.abs(ys - a*(xs - x0)**2 + y0))
    return x0, y0, 2*a, error # return min x, min y, second derivative


def minimize_starting_between(xlo, xhi, error_desired):
    xs = np.linspace(xlo, xhi, 5)
    es = np.array([func_to_minimize(x) for x in xs])
    plt.plot(xs,es,'.',label='data')
    x0, e0, deriv2, error = parabola_fit(xs, es)
    all_xs = np.linspace(xlo, xhi, 1000)
    print('x0', x0)
    print('e0', e0)
    plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
    plt.legend()
    plt.pause(1e-9)

    while True:
        # pick the next x:
        newx = x0 + np.random.normal()*np.sqrt(error/deriv2)
        xs = np.array(list(xs) + [newx])
        es = np.array(list(es) + [func_to_minimize(newx)])

        N_desired = max(5, (error/error_desired)**2)
        print('N_desired', N_desired)
        while len(xs) > N_desired:
            furthest = np.argmax(np.abs(xs-x0))
            xs[furthest] = xs[-1]
            es[furthest] = es[-1]
            xs = xs[:-1]
            es = es[:-1]

        plt.clf()
        plt.plot(xs,es,'.',label='data')
        x0, e0, deriv2, error = parabola_fit(xs, es)
        all_xs = np.linspace(xlo, xhi, 1000)
        print('x0', x0)
        print('e0', e0)
        plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
        plt.legend()
        plt.pause(1e-9)



minimize_starting_between(2.8, 3.6, 0.1)
