#!/usr/bin/python3

#Use this program to find the minimum of a set of data points subject to error with a Gaussian distribution
#Run this program from deft/papers/fuzzy-fmt with command python3 minimize_FE.py --kT [REQUIRED: temp]  --n [REQUIRED: density]   --tensor [optional]
#type --help for other options and default values

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
                    help='reduced density - REQUIRED')

#parser.add_argument('directory', metavar='directory', type=str,
                    #help='directory to store data - REQUIRED')

parser.add_argument('--numgw', default=10, type=float, help='number of inital gw datapoints')
parser.add_argument('--maxgw', type=float,
                    help='max gw', default=0.01)
parser.add_argument('--mingw', type=float,
                    help='min gw', default=0.001)

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
                    
parser.add_argument('--error_desired', metavar='error_desired', type=float,
                    help='error desired to be within', default=0.01)

parser.add_argument('--tensor', action='store_true',
                    help='use tensor weight')

args=parser.parse_args()

kT=args.kT
n=args.n


if args.maxgw:
    maxgw=args.maxgw
else :
    maxgw=0.6
    #maxgw=0.19   #needs to be much higher now!
    
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

def func_exact(x):
    return emin + 2*eps*((xmin/x)**12/2 - (xmin/x)**6) + eps
    
#def func_to_minimize(x):
    #return func_exact(x) + np.random.normal()*sigma

def func_to_minimize(kT, n, x, fv, dx, mcerror, mcconstant, mcprefactor):
    # return func_exact(x) + np.random.normal()*0.01
    print(kT,n,x,fv,dx,mcerror,mcconstant,mcprefactor)
    #sprintf(alldat_filedescriptor, "kT%5.3f_n%05.3f_fv%04.2f_gw%04.8f",
    name= 'kT%5.3f_n%05.3f_fv%04.2f_gw%04.8f' % (kT, n, fv, x)
    cmd = ' figs/new-melting.mkdat --kT %g --n %g' % (kT, n)
    cmd += ' --gw %.8g' % (x)
    cmd += ' --fv %g --dx %g' % (fv, dx)
    cmd += ' --mc-error %g --mc-constant %g --mc-prefactor %g' % (mcerror, mcconstant, mcprefactor)
    if args.tensor:
        os.system('mkdir -p newdata_tensor')
        cmd += ' --d newdata_tensor/phase-diagram'
        cmd += ' --filename isotherm-kT-%g-tensor.dat' % kT
    else:
        os.system('mkdir -p newdata')
        cmd += ' --d newdata/phase-diagram'
        cmd += ' --filename isotherm-kT-%g.dat' % kT
    if args.tensor:
        cmd += ' --tensor'
    print(cmd)
    os.system(cmd)
    if args.tensor:
        f='newdata_tensor/phase-diagram/%s-alldat_tensor.dat' % (name)
    else :
        f='newdata/phase-diagram/%s-alldat.dat' % (name)
    #if args.tensor:
        #f='%s/%s-alldat_tensor.dat' % (args.directory, name)
    #else :
        #f='%s/%s-alldat.dat' % (args.directory, name)
    data = np.loadtxt(f)   #FIX if gets a NAN, no data file created?
    diffFE=data[6]
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
    total_computations = args.numgw
    xs = np.linspace(xlo, xhi, args.numgw)  #steps by (xhi-xlo)/(total_computations-1)
    print('xs=',xs)
    #es = np.array([func_to_minimize(x) for x in xs])
    es = np.array([func_to_minimize(kT, n, x, fv, dx, mcerror, mcconstant, mcprefactor) for x in xs])
    x0, e0, deriv2, error = parabola_fit(xs, es)
    plt.plot(xs,es,'.',label='data')
    all_xs = np.linspace(xlo, xhi, 1000)
    print('gw', x0)
    print('diffFE', e0)
    plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
    lowest_indices = np.argsort(es)[:4]
    xs = xs[lowest_indices]
    es = es[lowest_indices]
    xs_to_fit = xs
    es_to_fit = es
    plt.plot(xs,es,'x',label='data to keep')
    plt.xlabel('gw')
    plt.ylabel('diffFE')
    plt.title("kT=%g, n=%g" % (kT, n))
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
        error_est = np.abs(rough_error_estimate(xs_to_fit,es_to_fit))
        if len(xs_to_fit) < 5:
            error_est=np.abs(rough_error_estimate(xs,es))

        newx = true_min_x + np.random.normal()*(xhi - xlo)/2
        xs = np.array(list(xs) + [newx])
        #es = np.array(list(es) + [func_to_minimize(newx)])
        es = np.array(list(es) + [func_to_minimize(kT, n, newx, fv, dx, mcerror, mcconstant, mcprefactor)])
        total_computations += 1

        xs_to_fit = []
        es_to_fit = []
        all_done = False
        while not all_done:
            best_height=5*error_est+height_y #matches minimum, but not curviture
            if np.isnan(best_height):
                print("crazy nan", error_est, height_y)
                exit(1)
            xlo = 1e300
            xhi = -1e300
            for i in range(len(xs)):
                if es[i] < best_height:
                    if xs[i] < xlo:
                        xlo = xs[i]
                    if xs[i] > xhi:
                        xhi = xs[i]

            for i in range(len(xs)):
                if xs[i] >= xlo and xs[i] <= xhi:
                    xs_to_fit.append(xs[i])
                    es_to_fit.append(es[i])

            if len(xs_to_fit) < 5:
                error_est *= 2
                xs_to_fit = []
                es_to_fit = []
            else:
                all_done = True
        xs_to_fit = np.array(xs_to_fit)
        es_to_fit = np.array(es_to_fit)

        x0, e0, deriv2, error = parabola_fit(xs_to_fit, es_to_fit)
        parabola_e = 0.5*deriv2*(xs - x0)**2 + e0
        residuals = np.abs(parabola_e - es)

        if xhi == true_min_x:
            xhi = xhi + 4*(xhi - xlo)
        if xlo == true_min_x:
            xlo = xlo - 4*(xhi - xlo)
        if xhi < x0 and deriv2 > 0:
            xhi = x0 + (x0 - xhi)
        if x0 < xlo and deriv2 > 0:
            xlo = x0 - (xlo - x0)

        num_samples = min(len([x for x in xs_to_fit if x < x0]),
                          len([x for x in xs_to_fit if x > x0]))
        error_estimate_of_min = 3*rough_error_estimate(xs_to_fit,es_to_fit)/np.sqrt(num_samples + 1e-10)

        plt.clf()
        plt.xlabel('gw')
        plt.ylabel('diffFE')
        if args.tensor:
            plt.title("kT=%g, n=%g, Tensor" % (kT, n))
        else:
            plt.title("kT=%g, n=%g, Nontensor" % (kT, n))
        plt.plot(xs,es,'.',label='data')
        plt.plot(xs_to_fit,es_to_fit,'.',label='data fitted', color='orange')
        all_xs = np.linspace(xs.min(), xs.max(), 1000)
        print('gw', x0, 'vs', true_min_x)
        print('diffFE', e0)
        plt.plot(all_xs, func_exact(all_xs), ':', label='exact')
        plt.plot(all_xs, 0.5*deriv2*(all_xs - x0)**2 + e0, label='my parabola')
        plt.axvline(true_min_x, color="blue")
        plt.axvline(xlo, color="red")
        plt.axvline(xhi, color="red")
        plt.axhline(e0 + error_estimate_of_min)
        plt.axhline(e0 - error_estimate_of_min)
        plt.axhline(best_height, linestyle=':')
        plt.legend()
        plt.xlim(min(xs_to_fit.min(), xlo) - 0.2, max(xs_to_fit.max(), xhi) + 0.2)
        plt.ylim(es_to_fit.min() - 0.2, es_to_fit.max() + 0.2)
        plt.pause(1.02)

        if error_estimate_of_min < error_desired and deriv2 > 0 and len(xs_to_fit) >= 5:
            print('excellent error:', error_estimate_of_min, 'with', len(xs), 'data and', total_computations, 'effort')
            print('rough_error_estimate', rough_error_estimate(xs_to_fit, es_to_fit), 'from', len(xs_to_fit))
            plt.show()
            return e0
        #plt.show()


minimize_starting_between(mingw, maxgw, args.error_desired)


