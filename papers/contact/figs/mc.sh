#!/bin/sh
#SBATCH --mem-per-cpu=2000
# #SBATCH --mail-type ALL
# #SBATCH --mail-user daveroundy@gmail.com
#SBATCH --output mc-%j.out

set -ev

hostname
date

MYNAME=`printf mc-%02d-%03d.dat $R $N`
MYTMP=`mktemp -d --tmpdir $MYNAME.XXXXXX`

echo time nice -19 ./monte-carlo $R $N $NITER `printf $MYTMP/mc-%02d-%03d.dat $R $N`
time nice -19 ./monte-carlo $R $N $NITER `printf $MYTMP/mc-%02d-%03d.dat $R $N`

ls -lh `printf $MYTMP/mc-%02d-%03d.dat $R $N`

echo time nice -19 ./analyze-monte-carlo $PREC `printf $MYTMP/mc-%02d-%03d.dat $R $N` `printf papers/contact/figs/mc-%02d-%03d.dat $R $N`
time nice -19 ./analyze-monte-carlo $PREC `printf $MYTMP/mc-%02d-%03d.dat $R $N` `printf papers/contact/figs/mc-%02d-%03d.dat $R $N`

rm `printf $MYTMP/mc-%02d-%03d.dat $R $N`
rmdir $MYTMP

date
