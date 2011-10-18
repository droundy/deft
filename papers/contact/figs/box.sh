#!/bin/sh
#SBATCH --mem-per-cpu=2000
# #SBATCH --mail-type ALL
# #SBATCH --mail-user daveroundy@gmail.com
#SBATCH --output box-%j.out

set -ev

hostname
date

if test -n "$XMAX"; then
    if test -n "$NUMBER"; then
        time nice -19 papers/contact/figs/box.mkdat $XMAX $YMAX $ZMAX $NUMBER
        date
        exit 0
    fi
fi

time nice -19 papers/contact/figs/sphere.mkdat
date
