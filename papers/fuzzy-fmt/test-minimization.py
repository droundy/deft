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

def func_number_datapoints(xs, xmin, xmax, es, ymin, ymax):
    count=0
    for i in range(1,len(xs)-1): 
        if xmin < xs[i] < xmax :
            if ymin < es[i] < ymax :
                count=count+1
    return count

def func_xcenter(xs, total_number_datapoints):
    sum_xs=0
    for i in range(1,len(xs)-1): 
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
        print("true_min", true_min)
        true_min_x = xs[true_min]
        print("true_min_x", true_min_x)
        height_y=es[true_min]
        print("height_y", height_y)
        best_height=5*np.abs(rough_error_estimate(xs,es))+height_y #matches minimum, but not curviture
        #best_height=10*np.abs(rough_error_estimate(xs,es))+height_y #good but shifted a little, closer curviture
        #best_height=20*np.abs(rough_error_estimate(xs,es))+height_y  
        print("best_height", best_height)
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
            #datapoints_on_left=func_number_datapoints(xs_to_fit,true_min_x-best_width, true_min_x, es_to_fit, height_y, best_height)
            datapoints_on_left=func_number_datapoints(xs,true_min_x-best_width, true_min_x, es, height_y, best_height)
            print("Number of datapoints on left side", datapoints_on_left)
            #datapoints_on_right=func_number_datapoints(xs_to_fit,true_min_x, true_min_x+best_width, es_to_fit, height_y, best_height)
            datapoints_on_right=func_number_datapoints(xs,true_min_x, true_min_x+best_width, es, height_y, best_height)
            print("Number of datapoints on right side", datapoints_on_right)
            #if datapoints_on_right > datapoints_on_left+10:  #NOT WORKING YET
                #true_min_x=true_min_x+func_xcenter(xs,datapoints_on_left+datapoints_on_right)
                #true_min_x=0.001*best_width+true_min_x
                #best_width_left=.8*best_width
                #print("best_width_left", best_width_left)
                #best_width_right=2*best_width-best_width_left
                #print("best_width_right", best_width_right)
            #if datapoints_on_left > datapoints_on_right+10:
                #true_min_x=true_min_x+func_xcenter(xs,datapoints_on_left+datapoints_on_right)
                #true_min_x=0.001*best_width-true_min_x
            furthest = np.argmax(np.abs(xs_to_fit-true_min_x))
            furthest_y = np.argmax(np.abs(es_to_fit-height_y))
            if np.abs(xs_to_fit[furthest] - true_min_x) > best_width:
                del xs_to_fit[furthest]
                del es_to_fit[furthest]
            elif np.abs(es_to_fit[furthest_y]) > best_height:
                del xs_to_fit[furthest_y]
                del es_to_fit[furthest_y]
            else :
                break

        xs_to_fit = np.array(xs_to_fit)
        es_to_fit = np.array(es_to_fit)
        x0, e0, deriv2, error = parabola_fit(xs_to_fit, es_to_fit)
        parabola_e = 0.5*deriv2*(xs - x0)**2 + e0
        residuals = np.abs(parabola_e - es)

        plt.clf()
        plt.plot(xs,es,'.',label='data')
        plt.plot(xs_to_fit,es_to_fit,'.',label='data fitted', color='orange')
        all_xs = np.linspace(xlo, xhi + (xhi-xlo)*.5, 1000)
        print('x0', x0, 'vs', true_min_x)
        print('e0', e0)
        plt.plot(all_xs, func_exact(all_xs), ':', label='exact')
        plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
        plt.axvline(true_min_x - best_width)
        plt.axvline(true_min_x + best_width)
        plt.axvline(true_min_x, color="blue")
        plt.axhline(y=height_y, color="blue")
        plt.axhline(y=best_height)
        plt.legend()
        plt.xlim(xs_to_fit.min() - 0.1, xs_to_fit.max() + 0.1)
        plt.ylim(es_to_fit.min() - 0.1, es_to_fit.max() + 0.1)
        #plt.pause(.2)
        plt.pause(.02)
        #plt.show()

        # xs = list(xs)
        # es = list(es)
        # for i in reversed(range(len(xs))):
        #     if residuals[i] > 3*error:
        #         del xs[i]
        #         del es[i]




minimize_starting_between(2.5, 3.8, 0.1)


