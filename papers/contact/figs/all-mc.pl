#!/usr/bin/perl -w

use strict;

my $iters  = 100000000;
my $acc = 0.001;
my $dir = "papers/contact/figs";

my $dd;
foreach $dd ([8,265], [6,112], [4,13]) {
    my $radius = sprintf("%02.0f", $$dd[0]);
    my $N = sprintf("%03.0f", $$dd[1]);

    print "Submitting monte-carlo in cavity with radius $radius\n";

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 0.001*($N) + 30; # It's a very hokey guess

    my $scriptname = "papers/contact/figs/mc-sphere-$radius-$N.tmp.sh";
    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
#SBATCH --mail-type FAIL
#SBATCH --mail-user daveroundy\@gmail.com
#SBATCH --output mc-sphere-$radius-$N.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 ./monte-carlo $N $iters $acc $dir/mc-$radius-$N.dat outerSphere $radius

date

";
    close(SCRIPT);
    system("sbatch", $scriptname);
}
