#!/bin/sh
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type ALL
#SBATCH --mail-user daveroundy@gmail.com
#SBATCH --output sphere-%j.out

set -ev

hostname
date

if test -n "$DIAMETER"; then
    if test -n "$NUMBER"; then
        time nice -19 papers/contact/figs/sphere.mkdat $DIAMETER $NUMBER
        date
        exit 0
    fi
fi

time nice -19 papers/contact/figs/sphere.mkdat
date
