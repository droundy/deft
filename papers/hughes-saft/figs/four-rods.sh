#!/bin/sh

for d in 0.6; do
    DIAMETER=$d sbatch papers/water-SAFT/figs/four-rods-in-water.sh
done
