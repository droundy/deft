#!/bin/sh
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type ALL
#SBATCH --mail-user jrbhughes@gmail.com

set -ev

hostname
date

time nice -19 paper/figs/single-rod-in-water.mkdat

date