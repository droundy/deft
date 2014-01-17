#!/usr/bin/perl -w

use strict;

my $t;
my $a;
# We do the largest diameters first, so the small calculations won't
# be scheduled on big-memory nodes leaving the big ones to wait for
# them to finish.
foreach $a ( "Kr" ) { #"Ne", "Ar", "Kr",  "Xe" ) {
  foreach $t ( 278 ) {
    my $atom = sprintf("%s", $a);
    my $temperature = sprintf("%04.2f", $t);

    print "Submitting $atom atom for $temperature K\n";
    my $padding = 1; # amount of extra space in nm that is added...
    my $resolution = 0.2; # spacing of grid points in nm.

    # Here I estimate the amount of memory that will be needed...
    my $memuse = 3000; #sprintf "%.0f", 2750*((($d + 2*$padding)/1.9)*(0.2/$resolution))**3;

    my $scriptname = "papers/water-saft/figs/hughes-lj-$atom-$temperature.tmp.sh";
    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user daveroundy\@gmail.com
#SBATCH --output hughes-lj-$a-$t.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 papers/water-saft/figs/hughes-lj-atom.mkdat $atom $temperature

date

";
  close(SCRIPT);
  system("sbatch", $scriptname);
}
}
