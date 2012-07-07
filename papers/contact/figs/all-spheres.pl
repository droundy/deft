#!/usr/bin/perl -w

use strict;

system "make papers/contact/figs/sphere.mkdat";

my $dd;
my $F;
# We do the largest radii first, so the small calculations won't
# be scheduled on big-memory nodes leaving the big ones to wait for
# them to finish.
foreach $dd ([20,0.5]#, [20,0.4], [20,0.3], [20,0.2], [20,0.1],
             #[16,0.5], [16,0.4], [16,0.3], [16,0.2], [16,0.1],
             #[12,0.5], [12,0.4], [12,0.3], [12,0.2], [12,0.1],
             #[ 8,0.5], [ 8,0.4], [ 8,0.3], [ 8,0.2], [ 8,0.1],
             #[ 6,0.5], [ 6,0.4], [ 6,0.3], [ 6,0.2], [ 6,0.1],
             #[ 4,0.5], [ 4,0.4], [ 4,0.3], [ 4,0.2], [ 4,0.1]
    ) {
  foreach $F ("WB", "WBm2") {
    my $r = sprintf("%02.0f", $$dd[0]);
    my $mu = sprintf("%04.1f", $$dd[1]);

    print "Submitting sphere in cavity $r\n";
    my $padding = 4; # amount of extra space in radii that is added...
    my $resolution = 0.1; # spacing of grid points in hard-sphere radius.

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 0.00025*((2*$r + $padding)/$resolution)**3;

    my $scriptname = "papers/contact/figs/outer-sphere-$r-$mu-$F.tmp.sh";
    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user daveroundy\@gmail.com
#SBATCH --output outer-sphere-$r-$mu-$F.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 papers/contact/figs/sphere.mkdat $r $mu $F

date

";
    close(SCRIPT);
    system("sbatch", $scriptname);
  }
}

