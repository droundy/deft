#!/bin/sh

# for d in 0.1 0.3 0.5 1.0; do
for d in `seq 0.1 0.1 0.9` 1.0 ; do
    #DIAMETER=$d ./sphere.sh
    DIAMETER=$d sbatch sphere.sh
done