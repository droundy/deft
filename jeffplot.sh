#!/bin/sh

#SBATCH --mem-per-cpu=1000


set -ev

python ./papers/contact/figs/plot-walls.py papers/contact/figs/mc-walls-20-19.dat papers/contact/figs/wallsWB-0.01.dat papers/contact/figs/wallsWBT-0.01.dat papers/contact/figs/wallsWBm2-0.01.dat jeffplot-0.01.pdf &

python ./papers/contact/figs/plot-walls.py papers/contact/figs/mc-walls-20-95.dat papers/contact/figs/wallsWB-0.05.dat papers/contact/figs/wallsWBT-0.05.dat papers/contact/figs/wallsWBm2-0.05.dat jeffplot-0.05.pdf &

python ./papers/contact/figs/plot-walls.py papers/contact/figs/mc-walls-20-193.dat papers/contact/figs/wallsWB-0.10.dat papers/contact/figs/wallsWBT-0.10.dat papers/contact/figs/wallsWBm2-0.10.dat jeffplot-0.10.pdf &

python ./papers/contact/figs/plot-walls.py papers/contact/figs/mc-walls-20-390.dat papers/contact/figs/wallsWB-0.20.dat papers/contact/figs/wallsWBT-0.20.dat papers/contact/figs/wallsWBm2-0.20.dat jeffplot-0.20.pdf &

python ./papers/contact/figs/plot-walls.py papers/contact/figs/mc-walls-20-589.dat papers/contact/figs/wallsWB-0.30.dat papers/contact/figs/wallsWBT-0.30.dat papers/contact/figs/wallsWBm2-0.30.dat jeffplot-0.30.pdf &

python ./papers/contact/figs/plot-walls.py papers/contact/figs/mc-walls-20-790.dat papers/contact/figs/wallsWB-0.40.dat papers/contact/figs/wallsWBT-0.40.dat papers/contact/figs/wallsWBm2-0.40.dat jeffplot-0.40.pdf &

