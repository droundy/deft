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

#def func_number_datapoints(xs, xmin, xmax, es, ymin, ymax):
def func_number_datapoints(xs, xmin, xmax):
    count=0
    for i in range(0,len(xs)):
        if xmin < xs[i] < xmax :
            #if ymin < es[i] < ymax :
            count=count+1
    return count

def func_xcenter(xs, xmin, xmax, total_number_datapoints):
    sum_xs=0
    for i in range(0,len(xs)):
        if  xmin < xs[i] < xmax :
            sum_xs=xs[i]+sum_xs
    return sum_xs/total_number_datapoints

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
    total_computations = 10
    xs = np.linspace(xlo, xhi, 10)
    es = np.array([func_to_minimize(x) for x in xs])
    x0, e0, deriv2, error = parabola_fit(xs, es)
    plt.plot(xs,es,'.',label='data')
    all_xs = np.linspace(xlo, xhi, 1000)
    print('x0', x0)
    print('e0', e0)
    plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
    lowest_indices = np.argsort(es)[:4]
    xs = xs[lowest_indices]
    es = es[lowest_indices]
    xs_to_fit = xs
    es_to_fit = es
    plt.plot(xs,es,'x',label='data to keep')
    plt.legend()
    plt.pause(0.4)

    while True:
        # pick the next x:
        true_min = np.argmin(np.abs(es))
        print("true_min", true_min, 'total effort', total_computations)
        true_min_x = xs[true_min]
        print("true_min_x", true_min_x)
        height_y=es[true_min]
        print("height_y", height_y)
        best_height=5*np.abs(rough_error_estimate(xs_to_fit,es_to_fit))+height_y #matches minimum, but not curviture
        if len(xs_to_fit) < 5:
            best_height=5*np.abs(rough_error_estimate(xs,es))+height_y
        print("best_height", best_height)
        newx = true_min_x + np.random.normal()*(xhi - xlo)/2
        xs = np.array(list(xs) + [newx])
        es = np.array(list(es) + [func_to_minimize(newx)])
        total_computations += 1

        N_desired = 10 # max(5, (error/error_desired)**2)
        print('N_desired', N_desired)
        N_desired = 1

        xlo = 1e300
        xhi = -1e300
        for i in range(len(xs)):
            if es[i] < best_height:
                if xs[i] < xlo:
                    xlo = xs[i]
                if xs[i] > xhi:
                    xhi = xs[i]
        if xhi == xlo:
            xhi = xs.max()
            xlo = xs.min()

        xs_to_fit = []
        es_to_fit = []
        for i in range(len(xs)):
            if xs[i] >= xlo and xs[i] <= xhi:
                xs_to_fit.append(xs[i])
                es_to_fit.append(es[i])
        xs_to_fit = np.array(xs_to_fit)
        es_to_fit = np.array(es_to_fit)

        x0, e0, deriv2, error = parabola_fit(xs_to_fit, es_to_fit)
        parabola_e = 0.5*deriv2*(xs - x0)**2 + e0
        residuals = np.abs(parabola_e - es)
        
        error_estimate_of_min = 3*rough_error_estimate(xs_to_fit,es_to_fit)/np.sqrt(len(xs))

        plt.clf()
        plt.plot(xs,es,'.',label='data')
        plt.plot(xs_to_fit,es_to_fit,'.',label='data fitted', color='orange')
        all_xs = np.linspace(xs.min(), xs.max(), 1000)
        print('x0', x0, 'vs', true_min_x)
        print('e0', e0)
        plt.plot(all_xs, func_exact(all_xs), ':', label='exact')
        plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
        plt.axvline(true_min_x, color="blue")
        plt.axhline(e0 + error_estimate_of_min)
        plt.axhline(e0 - error_estimate_of_min)
        plt.axhline(best_height, linestyle=':')
        plt.legend()
        plt.xlim(xs_to_fit.min() - 0.2, xs_to_fit.max() + 0.2)
        plt.ylim(es_to_fit.min() - 0.2, es_to_fit.max() + 0.2)
        plt.pause(1.02)

        if error_estimate_of_min < error_desired and deriv2 > 0:
            print('excellent error:', error_estimate_of_min, 'with', len(xs), 'data and', total_computations, 'effort')
            plt.show()
            return e0
        #plt.show()

        # xs = list(xs)
        # es = list(es)
        # for i in reversed(range(len(xs))):
        #     if residuals[i] > 3*error:
        #         del xs[i]
        #         del es[i]




minimize_starting_between(2.8, 3.8, 0.01)


