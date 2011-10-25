#!/bin/sh
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type ALL
#SBATCH --mail-user jrbhughes@gmail.com
#SBATCH --output single-rod-high-%j.out

set -ev

hostname
date

time nice -19 paper/figs/single-rod-in-water-high-res.mkdat

date