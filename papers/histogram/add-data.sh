#!/usr/bin/env sh

set -ev

min_seed=0
max_seed=29

min_N=5
max_golden_N=26
max_N=30

TMP=`mktemp -t tmp-XXXXXXXXX`

for N in `seq $min_N $max_N`; do
  for method in simple_flat tmmc oetmmc; do
    for ext in E lnw ps; do
      for seed in `seq $min_seed $max_seed`; do
        echo "data/s`printf '%03d' $seed`/periodic-ww1.30-ff0.30-N$N-$method-$ext.dat" >> $TMP
      done
    done
  done
done

for N in `seq $min_N $max_golden_N`; do
  for method in cfw wang_landau vanilla_wang_landau; do
    for ext in E lnw ps; do
      for seed in `seq $min_seed $max_seed`; do
        echo "data/s`printf '%03d' $seed`/periodic-ww1.30-ff0.30-N$N-$method-$ext.dat" >> $TMP
      done
    done
  done
done

for N in `seq $min_N $max_golden_N`; do
  for ext in E lnw ps transitions; do
    echo "data/periodic-ww1.30-ff0.30-N$N-tmmc-golden-$ext.dat" >> $TMP
  done
done

for seed in `seq $min_seed $max_seed`; do
  echo `ls data/s*$seed/periodic-ww1.30-ff0.30-N20-tmmc-transitions.dat` >> $TMP
done

git add -f `cat $TMP`

for seed in `seq $min_seed $max_seed`; do
  `dirname $0`/chseed.sh $seed;
done
`dirname $0`/chseed.sh 0
