#!/bin/sh

for d in 0.6 1.0 1.4 1.8 2.0 2.4; do
    DIAMETER=$d sbatch paper/figs/four-rods-in-water.sh
done
