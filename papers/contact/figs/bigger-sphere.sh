#!/bin/sh

set -ev

#!/bin/sh
#SBATCH --mem-per-cpu=3000
#SBATCH --mail-type ALL
#SBATCH --mail-user daveroundy@gmail.com
#SBATCH --output bigger-sphere-%j.out

set -ev

hostname
date

time nice -19 papers/contact/figs/sphere.mkdat 16 265 WB
date
