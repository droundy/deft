#!/usr/bin/perl -w

use strict;

system "make -j4 papers/water-saft/figs/hughes-single-rod.mkdat";

my $dd;
# We do the largest diameters first, so the small calculations won't
# be scheduled on big-memory nodes leaving the big ones to wait for
# them to finish.
#foreach $dd ( 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.9, 0.8,
#             0.7, 0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.17, 0.15, 0.13,
#             0.12, 0.11, 0.1, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0) {
foreach $dd ( 2.0, 1.6, 1.0, 0.6, 0.3, 0.1) {
  my $d = sprintf("%04.2f", $dd);

  print "Submitting single rod for $d nm\n";
  my $padding = 2; # amount of extra space in nm that is added...
  my $resolution = 0.1; # spacing of grid points in nm.

  # Here I estimate the amount of memory that will be needed...
  my $memuse = sprintf "%.0f", 425*((($d + 2*$padding)/5.0)*(0.1/$resolution))**2;

  my $scriptname = "papers/water-saft/figs/hughes-single-rod-$d.tmp.sh";
  open SCRIPT, ">$scriptname" or die $!;
  print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type END
##SBATCH --mail-user daveroundy\@gmail.com
#SBATCH --output single-$d.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 papers/water-saft/figs/hughes-single-rod.mkdat $d

date

";
  close(SCRIPT);
  #system("sh", $scriptname);
  system("sbatch", $scriptname);
}
