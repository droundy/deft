#!/bin/sh

set -ev

for N in `seq 8 2 14`; do
    minD=`echo "printf(\"%g\", ceil(1 + ($N/0.3)**(1/3)))" | octave -q`
    echo minD is $minD for N $N
    for D in `seq $minD 6`; do
        echo XMAX=$D YMAX=$D ZMAX=$D NUMBER=$N sbatch -J box-$D-$N papers/contact/figs/box.sh
        XMAX=$D YMAX=$D ZMAX=$D NUMBER=$N sbatch -J box-$D-$N papers/contact/figs/box.sh
    done
done
