#!/usr/bin/env sh

#set -ev

ff=0.3
ww=1.3
min_T=0.1

for N in `seq 5 30`; do
    echo Working on N=$N
    python2 run-golden.py $ww $ff $min_T $N
done
