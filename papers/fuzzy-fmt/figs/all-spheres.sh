#!/bin/sh

set -ev

DIAMETER=16 NUMBER=265 FUNCTIONAL=WB sbatch -J sphereWB-16-265 papers/contact/figs/sphere.sh
DIAMETER=16 NUMBER=265 FUNCTIONAL=WBT sbatch -J sphereWBT-16-265 papers/contact/figs/sphere.sh

DIAMETER=12 NUMBER=112 FUNCTIONAL=WB sbatch -J sphereWB-12-112 papers/contact/figs/sphere.sh
DIAMETER=12 NUMBER=112 FUNCTIONAL=WBT sbatch -J sphereWBT-12-112 papers/contact/figs/sphere.sh

for N in `seq 2 13`; do
    minD=`echo "printf(\"%g\", 2*ceil(1 + (3*$N/4/pi/0.5)**(1/3)))" | octave -q`
    echo minD is $minD for N $N
    for D in `seq $minD 8`; do
        for F in WB WBT; do
            echo DIAMETER=$D NUMBER=$N FUNCTIONAL=$F sbatch papers/contact/figs/sphere.sh
            DIAMETER=$D NUMBER=$N FUNCTIONAL=$F sbatch -J sphere$F-$D-$N papers/contact/figs/sphere.sh
        done
    done
done
