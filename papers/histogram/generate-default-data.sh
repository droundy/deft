#!/usr/bin/env sh

#set -ev

ff=0.3
ww=1.3
min_T=0.2

echo using method: $method
for N in `seq 5 30`; do
  echo '  working on N='$N
  for method in simple_flat wang_landau tmmc oetmmc; do
    echo working on $method
    python2 run-default.py $ww $ff $min_T $N $method
  done
done
