#!/usr/bin/perl -w

use strict;

system "make papers/contact/figs/sphere.mkdat";

my $dd;
my $F;
# We do the largest diameters first, so the small calculations won't
# be scheduled on big-memory nodes leaving the big ones to wait for
# them to finish.
foreach $dd ([20,4000], [20,2000], [16,2050], [16,850], [12,400], [12,860]) {
  foreach $F ("WB", "WBT", "WBm2") {
    my $d = sprintf("%02.0f", $$dd[0]);
    my $N = sprintf("%04.0f", $$dd[1]);

    print "Submitting sphere in cavity $d\n";
    my $padding = 4; # amount of extra space in radii that is added...
    my $resolution = 0.1; # spacing of grid points in hard-sphere radius.

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 2400*((($d + $padding)/20.0)*(0.1/$resolution))**3;

    my $scriptname = "papers/contact/figs/sphere-$d-$N-$F.tmp.sh";
    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user daveroundy\@gmail.com
#SBATCH --output sphere-$d-$N-$F.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 papers/contact/figs/sphere.mkdat $d $N $F

date

";
    close(SCRIPT);
    system("sbatch", $scriptname);
  }
}

