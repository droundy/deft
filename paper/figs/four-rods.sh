#!/bin/sh

for d in 2.0; do
    DIAMETER=$d sbatch paper/figs/four-rods-in-water.sh
done
