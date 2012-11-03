#!/bin/sh

for d in 0.4 0.8 1.2 1.6 2.0 2.4; do
    DIAMETER=$d papers/water-SAFT/figs/rods-in-water.sh
    #DIAMETER=$d sbatch papers/water-SAFT/figs/rods-in-water.sh
done
