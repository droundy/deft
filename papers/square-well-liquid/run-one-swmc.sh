#!/usr/bin/sh

# this script runs one method with default arguments and updates figures with the new data

set -ev

ww=1.3
ff=0.3
N=20
iterations=1000000

if [ $# -eq 1 ]; then
  cores=4
elif [ $# -eq 2 ]; then
  cores=$2
else
  echo "usage: [script] method cores"
  exit 1
fi

error_hook(){
  if ! [ $? -eq 0 ]; then
    echo "error with exit status" $?
    exit $?
  fi
}

method=$1

# move to the deft directory
cd `dirname $0`/../../
scons -j$cores square-well-monte-carlo
error_hook
./square-well-monte-carlo --ww $ww --ff $ff --N $N --iterations $iterations --$method
error_hook
./papers/square-well-liquid/build-swmc-paper.py

