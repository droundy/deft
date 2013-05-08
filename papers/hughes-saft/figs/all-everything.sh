#!/bin/sh

echo Submitting all slurm jobs needed to regenerate paper!

set -ev

papers/water-SAFT/figs/all-spheres.pl

papers/water-SAFT/figs/all-single-rods.pl

papers/water-SAFT/figs/all-rods.sh

papers/water-SAFT/figs/all-four-rods.sh

sbatch papers/water-SAFT/figs/single-rod-high.sh

sbatch papers/water-SAFT/figs/single-rod-low.sh
