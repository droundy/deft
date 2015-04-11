#!/usr/bin/env sh

#set -ev

ff=0.3
ww=1.3
min_T=0.2

for N in `seq 5 30`; do
  for seed in `seq 0 4`; do
    echo '  working on N='$N
    for method in simple_flat tmmc oetmmc; do
      #for method in cfw wang_landau; do
      echo working on $method
      python2 run-default.py $ww $ff $min_T $N $method $seed
    done
  done
done
