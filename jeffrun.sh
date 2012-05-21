#!/bin/sh

#SBATCH --mem-per-cpu=1000


set -ev

srun -J mc-periodicxyz-20-193 ./monte-carlo 193 300000000 .001 ./papers/contact/figs/mc-periodicxyz-20-193.dat periodxyz 20 flatdiv > mc-periodicxyz-20-193.out & 

srun -J mc-periodicxyz-20-390 ./monte-carlo 390 300000000 .001 ./papers/contact/figs/mc-periodicxyz-20-390.dat periodxyz 20 flatdiv > mc-periodicxyz-20-390.out &

srun -J mc-periodicxyz-20-589 ./monte-carlo 589 300000000 .001 ./papers/contact/figs/mc-periodicxyz-20-589.dat periodxyz 20 flatdiv > mc-periodicxyz-20-589.out &

srun -J mc-periodicxyz-20-790 ./monte-carlo 790 300000000 .001 ./papers/contact/figs/mc-periodicxyz-20-790.dat periodxyz 20 flatdiv > mc-periodicxyz-20-790.out &

srun -J mc-periodicxyz-20-95 ./monte-carlo 95 3000000000 .001 ./papers/contact/figs/mc-periodicxyz-20-95.dat periodxyz 20 flatdiv > mc-periodicxyz-20-95.out &

srun -J mc-periodicxyz-20-19 ./monte-carlo 19 3000000000 .001 ./papers/contact/figs/mc-periodicxyz-20-19.dat periodxyz 20 flatdiv > mc-periodicxyz-20-19.out &

srun -J mc-periodicxyz-20-134 ./monte-carlo 134 3000000000 .001 ./papers/contact/figs/mc-periodicxyz-20-134.dat periodxyz 20 flatdiv > mc-periodicxyz-20-134.out &

srun -J mc-walls-20-193 ./monte-carlo 193 300000000 .001 ./papers/contact/figs/mc-walls-20-193.dat periodxy 20 wallz 20 flatdiv > mc-walls-20-193.out &

srun -J mc-walls-20-95 ./monte-carlo 95 3000000000 .001 ./papers/contact/figs/mc-walls-20-95.dat periodxy 20 wallz 20 flatdiv > mc-walls-20-95.out &

srun -J mc-walls-20-19 ./monte-carlo 19 3000000000 .001 ./papers/contact/figs/mc-walls-20-19.dat periodxy 20 wallz 20 flatdiv > mc-walls-20-19.out &

srun -J mc-walls-20-390 ./monte-carlo 390 3000000000 .001 ./papers/contact/figs/mc-walls-20-390.dat periodxy 20 wallz 20 flatdiv > mc-walls-20-390.out &

srun -J mc-walls-20-589 ./monte-carlo 589 3000000000 .001 ./papers/contact/figs/mc-walls-20-589.dat periodxy 20 wallz 20 flatdiv > mc-walls-20-589.out &

srun -J mc-walls-20-790 ./monte-carlo 790 3000000000 .001 ./papers/contact/figs/mc-walls-20-790.dat periodxy 20 wallz 20 flatdiv > mc-walls-20-790.out &


srun -J mc-walls-23-146 ./monte-carlo 146 300000000 .001 ./papers/contact/figs/mc-walls-23-146.dat periodxy 20 wallz 20 flatdiv > mc-walls-23-146.out &

srun -J mc-walls-39-142 ./monte-carlo 142 3000000000 .001 ./papers/contact/figs/mc-walls-39-142.dat periodxy 20 wallz 20 flatdiv > mc-walls-39-142.out &

srun -J mc-walls-44-143 ./monte-carlo 44 3000000000 .001 ./papers/contact/figs/mc-walls-44-143.dat periodxy 20 wallz 20 flatdiv > mc-walls-44-143.out &

srun -J mc-periodicxyz-23-146 ./monte-carlo 146 3000000000 .001 ./papers/contact/figs/mc-periodicxyz-23-146.dat periodxyz 20 flatdiv > mc-periodicxyz-23-146.out &

srun -J mc-periodicxyz-39-142 ./monte-carlo 142 3000000000 .001 ./papers/contact/figs/mc-periodicxyz-39-142.dat periodxyz 20 flatdiv > mc-periodicxyz-39-142.out &

srun -J mc-periodicxyz-44-143 ./monte-carlo 143 3000000000 .001 ./papers/contact/figs/mc-periodicxyz-44-143.dat periodxyz 20 flatdiv > mc-periodicxyz-44-143.out &