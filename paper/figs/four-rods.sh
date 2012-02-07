#!/bin/sh

for d in 0.6; do
    DIAMETER=$d sbatch paper/figs/four-rods-in-water.sh
done
