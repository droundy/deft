#!/bin/sh
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type ALL
#SBATCH --mail-user daveroundy@gmail.com
#SBATCH --output sphere-%j.out

set -ev

hostname
date

if test -n "$DIAMETER"; then
    time nice -19 paper/figs/sphere.mkdat $DIAMETER
else
    time nice -19 paper/figs/sphere.mkdat
fi


date
