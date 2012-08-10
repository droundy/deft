#!/usr/bin/perl -w

use strict;

system "make papers/contact/figs/inner-sphere.mkdat";

my $dd;
my $F;
# We do the largest diameters first, so the small calculations won't
# be scheduled on big-memory nodes leaving the big ones to wait for
# them to finish.
foreach $dd ([4,0.1], [4,0.2], [4,0.3], [4,0.4], [4,0.5],
             [12,0.4], [12, 0.1]
            ) {
  foreach $F ("WB", "WBm2") {
    my $d = sprintf("%02.0f", $$dd[0]);
    my $eta = sprintf("%04.2f", $$dd[1]);

    print "Submitting sphere with diameter $d\n";
    my $padding = 16; # amount of extra space in radii that is added...
    my $resolution = 0.1; # spacing of grid points in hard-sphere radius.

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 2500*((($d + $padding)/20.0)*(0.1/$resolution))**3;

    my $scriptname = "papers/contact/figs/inner-sphere-$d-$eta-$F.tmp.sh";
    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user daveroundy\@gmail.com
#SBATCH --output inner-sphere-$d-$eta-$F.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 papers/contact/figs/inner-sphere.mkdat $d $eta $F

date

";
    close(SCRIPT);
    system("sbatch", $scriptname);
  }
}

