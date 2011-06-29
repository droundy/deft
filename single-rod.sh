#!/bin/sh
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type ALL
#SBATCH --mail-user daveroundy@gmail.com
#SBATCH --output single-rod-%j.out

set -ev

hostname
date

time nice -19 paper/figs/single-rod-in-water.mkdat

date
