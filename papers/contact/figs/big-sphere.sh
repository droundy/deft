#!/bin/sh

set -ev

#!/bin/sh
#SBATCH --mem-per-cpu=2000
#SBATCH --mail-type ALL
#SBATCH --mail-user daveroundy@gmail.com
#SBATCH --output big-sphere-%j.out

set -ev

hostname
date

time nice -19 papers/contact/figs/sphere.mkdat 12 112
date
