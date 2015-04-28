#!/usr/bin/env sh

#set -ev

ff=0.3
ww=1.3
min_T=0.2

min_N=5
max_N=30
max_golden_N=26

for N in `seq $min_N $max_N`; do
  echo "\n\n\n    working on N=$N\n"
  for seed in `seq 18 24`; do
    echo "\n\n  working on seed=$seed\n"

    for method in simple_flat tmmc oetmmc; do
      echo "\nworking on $method\n"
      python2 run-default.py $ww $ff $min_T $N $method $seed
    done

    if [ $N -le $max_golden_N ]; then
      for method in cfw wang_landau; do
        echo "\nworking on $method\n"
        python2 run-default.py $ww $ff $min_T $N $method $seed
      done
    fi

  done
done
