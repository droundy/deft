#!/bin/sh
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type ALL
#SBATCH --mail-user jrbhughes@gmail.com
#SBATCH --output rods-%j.out

set -ev

hostname
date

if test -n "$DIAMETER"; then
    time nice -19 paper/figs/rods-in-water.mkdat $DIAMETER
else
    time nice -19 paper/figs/rods-in-water.mkdat
fi

date