#!/bin/bash

set -ev

MC=papers/contact/figs/mc.sh

time make -j4 monte-carlo analyze-monte-carlo

rm -f mc-*.out

export NITER=10000000
export PREC=0.1

export R=8
export N=265
sbatch -J mc-$R-$N -o mc-$R-$N.out $MC

export R=6
export N=112
sbatch -J mc-$R-$N -o mc-$R-$N.out $MC

export R=4
export N=13
sbatch -J mc-$R-$N -o mc-$R-$N.out $MC
