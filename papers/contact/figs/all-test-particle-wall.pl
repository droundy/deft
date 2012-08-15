#!/usr/bin/perl -w

use strict;

system "make papers/contact/figs/test-particle-wall.mkdat";

my $dd;
my $F;
# We do the largest diameters first, so the small calculations won't
# be scheduled on big-memory nodes leaving the big ones to wait for
# them to finish.
foreach $dd ([2,0.4,3], [2,0.4,4], [2,0.4,5], [2,0.4,6], [2,0.4,8], [2,0.4,11],
    [2,0.1,3], [2,0.1,4], [2,0.1,5], [2,0.1,11]) {
  foreach $F ("WB"){
#, "WBT", "WBm2") {
    my $d = sprintf("%02.0f", $$dd[0]);
    my $eta = sprintf("%04.2f", $$dd[1]);
    my $z = sprintf("%04.2f", $$dd[2]);

    print "Submitting sphere with diameter $d that is $z away from origin\n";
    # Here I estimate the amount of memory that will be needed...
    my $memuse = 4000;

    my $scriptname = "papers/contact/figs/test-particle-wall-$d-$eta-$z-$F.tmp.sh";
    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user daveroundy\@gmail.com
#SBATCH --output test-particle-wall-$d-$eta-$z-$F.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 papers/contact/figs/test-particle-wall.mkdat $d $eta $z $F

date

";
    close(SCRIPT);
    system("sbatch", $scriptname);
  }
}
