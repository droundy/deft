#!/bin/sh

echo Submitting all slurm jobs needed to regenerate paper!

set -ev

paper/figs/all-spheres.pl

paper/figs/all-single-rods.pl

paper/figs/all-rods.sh

paper/figs/four-rods.sh

sbatch paper/figs/single-rod-high.sh

sbatch paper/figs/single-rod-low.sh
