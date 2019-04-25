#!/usr/bin/python3

#Use this program to find the minimum of a set of data points subject to error with a Gaussian distribution
#Run this program from deft/papers/fuzzy-fmt with command python3 test-minimization.py

from __future__ import division, print_function

import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

parser = argparse.ArgumentParser(description='Finds the minimum difference in FE.')

parser.add_argument('--kT', metavar='temperature', type=float,
                    help='reduced temperature - REQUIRED')
parser.add_argument('--n', type=float,
                    help='density - REQUIRED')

parser.add_argument('--numgw', type=float,
                    help='number of inital gw', default=10)
parser.add_argument('--maxgw', type=float,
                    help='max gw', default=0.2)
parser.add_argument('--mingw', type=float,
                    help='min gw', default=0.01)

parser.add_argument('--fv', metavar='vacancies', type=float,
                    help='fraction of vacancies - Default 0')
parser.add_argument('--dx', metavar='dx', type=float,
                    help='scaling dx - Default 0.5')
parser.add_argument('--mcerror', metavar='mc_error', type=float,
                    help='monte carlo mc_error - Default 0.001')
parser.add_argument('--mcconstant', metavar='const', type=int,
                    help='monte carlo integration mc_constant - Default 5')
parser.add_argument('--mcprefactor', metavar='prefac', type=int,
                    help='monte carlo integration mc_prefactor - Default 50000')
                    
parser.add_argument('--tensor', action='store_true',
                    help='use tensor weight')

args=parser.parse_args()

kT=args.kT
n=args.n


if args.numgw:
    numgw=args.numgw
else :
    dgw=10

if args.maxgw:
    maxgw=args.maxgw
else :
    maxgw=0.2  
    
if args.mingw:
    mingw=args.mingw
else :
    mingw=0.01

if args.fv:
    fv=args.fv
else :
    fv=0

if args.dx:
    dx=args.dx
else :
    dx=.5

if args.mcerror:
    mcerror=args.mcerror
else :
    mcerror=0.001

if args.mcconstant:
    mcconstant=args.mcconstant
else :
    mcconstant=5
    
if args.mcprefactor:
    mcprefactor=args.mcprefactor
else :
    mcprefactor=50000


xmin = 3.0
eps = 2.0
emin = 1.0
sigma = 0.01

#def func_number_datapoints(xs, xmin, xmax, es, ymin, ymax):
#def func_number_datapoints(xs, xmin, xmax):
    #count=0
    #for i in range(0,len(xs)):
        #if xmin < xs[i] < xmax :
            ##if ymin < es[i] < ymax :
            #count=count+1
    #return count

#def func_xcenter(xs, xmin, xmax, total_number_datapoints):
    #sum_xs=0
    #for i in range(0,len(xs)):
        #if  xmin < xs[i] < xmax :
            #sum_xs=xs[i]+sum_xs
    #return sum_xs/total_number_datapoints

def func_exact(x):
    return emin + 2*eps*((xmin/x)**12/2 - (xmin/x)**6) + eps
    
#def func_to_minimize(x):
    #return func_exact(x) + np.random.normal()*sigma

def func_to_minimize(kT, n, x, fv, dx, mcerror, mcconstant, mcprefactor):
    print(kT,n,x,fv,dx,mcerror,mcconstant,mcprefactor)
    cmd = ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
    cmd += ' --gw %g' % (x)
    cmd += ' --fv %g --dx %g' % (fv, dx)
    cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
    cmd += ' --filename isotherm-kT-%g.dat' % kT
    print(cmd)
    os.system(cmd)
    f='crystallization/kT%.3f_n%.3f_fv%.2f_gw%.3f-alldat.dat' % (kT, n, fv, x)
    data = np.loadtxt(f)
    diffFE=data[6]
    print('diffFE=', diffFE, 'gw=', x)
    print()
    print()
    return diffFE
    

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
    xs = np.linspace(xlo, xhi, 10)  #steps by (xhi-xlo)/(total_computations-1)
    print('xs=',xs)
    #es = np.array([func_to_minimize(x) for x in xs])
    es = np.array([func_to_minimize(kT, n, x, fv, dx, mcerror, mcconstant, mcprefactor) for x in xs])
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
        #es = np.array(list(es) + [func_to_minimize(newx)])
        es = np.array(list(es) + [func_to_minimize(kT, n, newx, fv, dx, mcerror, mcconstant, mcprefactor)])
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




minimize_starting_between(0.01, 0.19, 0.01)


