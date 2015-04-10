#!/usr/bin/env sh

#set -ev

ff=0.3
ww=1.3
min_T=0.2

for N in 70 80 90 100; do
  echo '  working on N='$N
  #python2 run-default.py $ww $ff $min_T $N transitions_input_filename data/periodic-ww1.30-ff0.30-N$N-tmmc-golden-transitions.dat
  for method in simple_flat tmmc oetmmc; do
    echo working on $method
    python2 run-default.py $ww $ff $min_T $N $method
  done
done
