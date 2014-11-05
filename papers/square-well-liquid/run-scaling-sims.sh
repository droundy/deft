#!/bin/sh

ww=1.3
ff=0.3

for N in 10 20 40 100 200; do
  python2 `dirname $0`/run-monte-carlo.py $ww $ff $N '["wang_landau","flat","gaussian"]' "scaling"
done
