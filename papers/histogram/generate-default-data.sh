#!/usr/bin/env sh

#set -ev

ff=0.3
ww=1.3
min_T=0.2

for N in `seq 27 30`; do
  echo "\n\n\n    working on N=$N\n"
  for seed in `seq 0 9`; do
    echo "\n\n  working on seed=$seed\n"
    #for method in simple_flat tmmc oetmmc; do
    for method in cfw wang_landau; do
      echo "\nworking on $method\n"
      python2 run-default.py $ww $ff $min_T $N $method $seed
    done
  done
done
