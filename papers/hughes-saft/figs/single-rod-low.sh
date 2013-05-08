#!/bin/sh
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type ALL
#SBATCH --mail-user jrbhughes@gmail.com
#SBATCH --output single-rod-low-%j.out

set -ev

hostname
date

time nice -19 papers/water-SAFT/figs/single-rod-in-water-low-res.mkdat

date
