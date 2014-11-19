#!/bin/sh

ww=1.3
ff=0.3
N=60

python2 `dirname $0`/run-monte-carlo.py $ww $ff $N '["wang_landau","flat","gaussian"]' "scaling"
