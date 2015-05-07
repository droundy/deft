#!/usr/bin/env sh

to=$1
TMP=`mktemp -t tmp-XXXXXXXXX`

for file in figs/*.py; do
  cat $file | sed "s/seed\ =\ \[.*\]/seed\ =\ \[$to\]/" > $TMP
  mv $TMP $file
done

git add -f `cat $TMP`
