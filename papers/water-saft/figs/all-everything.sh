#!/bin/sh

echo Submitting all slurm jobs needed to regenerate paper!

set -ev

papers/water-saft/figs/all-spheres.pl

papers/water-saft/figs/all-single-rods.pl

papers/water-saft/figs/all-rods.sh

papers/water-saft/figs/all-four-rods.sh

sbatch papers/water-saft/figs/single-rod-high.sh

sbatch papers/water-saft/figs/single-rod-low.sh
