#!/bin/sh

for d in 1.0; do
    #DIAMETER=$d ./rods-in-water.sh
    DIAMETER=$d sbatch rods-in-water.sh
done