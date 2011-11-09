#!/usr/bin/perl -w

use strict;

my $d;
# We do the largest diameters first, so the small calculations won't
# be scheduled on big-memory nodes leaving the big ones to wait for
# them to finish.
foreach $d ('1.0', 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, '0.0') {

  print "Submitting sphere for $d nm\n";
  my $padding = 1; # amount of extra space in nm that is added...
  my $resolution = 0.2; # spacing of grid points in nm.

  # Here I estimate the amount of memory that will be needed.  I know
  # from experience that 2000 works fine for a 1.9 nm cubical cell at
  # a resolution of 0.2 nm.
  my $memuse = sprintf "%.0f", 2000*((($d + 2*$padding)/1.9)*(0.2/$resolution))**3;

  my $scriptname = "paper/figs/sphere-$d.tmp.sh";
  open SCRIPT, ">$scriptname" or die $!;
  print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
#SBATCH --mail-type ALL
#SBATCH --mail-user daveroundy\@gmail.com
#SBATCH --output sphere-$d.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 paper/figs/sphere.mkdat $d

date

";
  close(SCRIPT);
  system("sbatch", $scriptname);
}
