#!/bin/sh

ww=1.3
ff=0.3

for N in 20 60 100 200; do
  python2 `dirname $0`/run-monte-carlo.py $ww $ff $N '["wang_landau","flat","gaussian"]' "scaling"
done
