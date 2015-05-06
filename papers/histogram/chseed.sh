#!/usr/bin/env sh

to=$1
TMP=`mktemp -t tmp-XXXXXXXXX`

for file in figs/*.py; do
  cat $file | sed "s/seed\ =\ \[.*\]/seed\ =\ \[$to\]/" > $TMP
  mv $TMP $file
done

cd `dirname $0`/data/s*$1

for ext in E lnw; do
  file="periodic-ww1.30-ff0.30-N20-tmmc-golden-$ext.dat"
  ln -s ../$file ./ 2>/dev/null
  echo $file >> $TMP
done

for ext in E lnw ps; do
  file="periodic-ww1.30-ff0.30-N20-nw-$ext.dat"
  ln -s ../s000/$file ./ 2>/dev/null
  echo $file >> $TMP
done

for kT in 0.4 0.5; do
  file="periodic-ww1.30-ff0.30-N20-kT$kT-E.dat"
  ln -s ../s000/$file ./ 2>/dev/null
  echo $file >> $TMP
done

git add -f `cat $TMP`
