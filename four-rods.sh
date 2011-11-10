#!/bin/sh

for d in 1.0; do
    DIAMETER=$d sbatch four-rods-in-water.sh
done
