#!/bin/sh

set -ev

for N in `seq 5 20`; do
    echo Working on N=$N
    nice -19 ./square-well-monte-carlo --N $N --filename_suffix "golden" --min_T 0.1 --min_samples 10000 --tmmc --ff 0.3 --ww 1.3 --iterations 100000000
done
