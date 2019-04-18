#!/usr/bin/python3

#Use this program to find the minimum of a set of data points 
#Run this program from deft/papers/fuzzy-fmt with command python3 test-minimization.py

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

xmin = 3.0
eps = 2.0
emin = 1.0
sigma = 0.01

def func_exact(x):
    return emin + 2*eps*((xmin/x)**12/2 - (xmin/x)**6) + eps

def func_to_minimize(x):
    return func_exact(x) + np.random.normal()*sigma

def rough_error_estimate(xs,ys):
    iii = np.argsort(xs)
    xs = np.array(xs[iii])
    ys = np.array(ys[iii])
    total = 0
    for i in range(1,len(xs)-1):
        total += abs(ys[i] - ((xs[i] - xs[i-1])*ys[i+1] + (xs[i+1]-xs[i])*ys[i-1])/(xs[i+1]-xs[i-1]))
    return total/(len(xs)-2)

def parabola_fit(xs,ys):
    A = np.vstack([xs**2, xs, np.ones(len(xs))]).T
    a, b, c = np.linalg.lstsq(A, ys, rcond=-1)[0]
    x0 = -b/(2*a)
    y0 = a*x0**2 + b*x0 + c
    error = np.mean(np.abs(ys - a*(xs - x0)**2 + y0))
    print('error', error)
    return x0, y0, 2*a, error  # return min x, min y, second derivative


def minimize_starting_between(xlo, xhi, error_desired):
    xs = np.linspace(xlo, xhi, 5)
    es = np.array([func_to_minimize(x) for x in xs])
    x0, e0, deriv2, error = parabola_fit(xs, es)
    plt.plot(xs,es,'.',label='data')
    all_xs = np.linspace(xlo, xhi, 1000)
    print('x0', x0)
    print('e0', e0)
    plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
    plt.legend()
    plt.pause(1e-9)
    best_width = np.sqrt(error/deriv2)

    while True:
        # pick the next x:
        true_min = np.argmin(np.abs(es))
        true_min_x = xs[true_min]
        if deriv2 > 0:
            best_width = min(np.sqrt(error/deriv2), best_width)
        newx = np.random.random()*(2*best_width) + true_min_x - best_width
        xs = np.array(list(xs) + [newx])
        es = np.array(list(es) + [func_to_minimize(newx)])

        N_desired = 10 # max(5, (error/error_desired)**2)
        print('N_desired', N_desired)
        N_desired = 1


        xs_to_fit = list(1*xs)
        es_to_fit = list(1*es)
        while len(xs_to_fit) > N_desired:
            furthest = np.argmax(np.abs(xs_to_fit-true_min_x))
            if np.abs(xs_to_fit[furthest] - true_min_x) > best_width:
                del xs_to_fit[furthest]
                del es_to_fit[furthest]
            else:
                break
        xs_to_fit = np.array(xs_to_fit)
        es_to_fit = np.array(es_to_fit)
        x0, e0, deriv2, error = parabola_fit(xs_to_fit, es_to_fit)
        parabola_e = 0.5*deriv2*(xs - x0)**2 + e0
        residuals = np.abs(parabola_e - es)

        plt.clf()
        plt.plot(xs,es,'.',label='data')
        plt.plot(xs_to_fit,es_to_fit,'.',label='data fitted')
        all_xs = np.linspace(xlo, xhi + (xhi-xlo)*.5, 1000)
        print('x0', x0, 'vs', true_min_x)
        print('e0', e0)
        plt.plot(all_xs, func_exact(all_xs), ':', label='exact')
        plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
        plt.axvline(true_min_x - best_width)
        plt.axvline(true_min_x + best_width)
        plt.legend()
        plt.xlim(xs_to_fit.min() - 0.1, xs_to_fit.max() + 0.1)
        plt.ylim(es_to_fit.min() - 0.1, es_to_fit.max() + 0.1)
        plt.pause(1)

        # xs = list(xs)
        # es = list(es)
        # for i in reversed(range(len(xs))):
        #     if residuals[i] > 3*error:
        #         del xs[i]
        #         del es[i]




minimize_starting_between(2.5, 3.8, 0.1)
