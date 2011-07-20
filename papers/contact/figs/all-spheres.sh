#!/bin/sh

set -ev

for N in `seq 2 13`; do
    minD=`echo "printf(\"%g\", 2*ceil(1 + (3*$N/4/pi/0.5)**(1/3)))" | octave -q`
    echo minD is $minD for N $N
    for D in `seq $minD 8`; do
        echo DIAMETER=$D NUMBER=$N sbatch papers/contact/figs/sphere.sh
        DIAMETER=$D NUMBER=$N sbatch -J sphere-$D-$N papers/contact/figs/sphere.sh
    done
done
