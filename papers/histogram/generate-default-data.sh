#!/usr/bin/env sh

#set -ev

ff=0.3
ww=1.3
min_T=0.2
method=simple_flat
# method is one of: simple_flat wang_landau tmmc oetmmc

echo using method: $method
for N in `seq 5 20` 40 60 80 100 120 140 160 180 200 250 300; do
  echo working on N=$N
  python2 run-default.py $ww $ff $min_T $N $method
done
