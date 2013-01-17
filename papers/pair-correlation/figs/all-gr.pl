#!/usr/bin/perl -w

use strict;
use Digest::SHA1 qw(sha1_hex);
use File::Slurp qw(read_file);

my $dir = "papers/pair-correlation/figs";
my $sha1code = sha1_hex(read_file("src/Monte-Carlo/radial-distribution-monte-carlo.cpp"));

system "make radial-distribution-monte-carlo" || die "make failed";

my $acc = 1e-5; # This is the fractional uncertainty in g(r) at contact
my $dr = 0.01;
my $version = "";

my $dd;
foreach $dd (0.1, 0.2, 0.3, 0.4, 0.5) {
    my $packingfraction = sprintf("%04.2f", $dd);

    # Here I estimate the amount of memory that will be needed...
    my $memuse = 30; # It's a very hokey guess

    my $scriptname = "$dir/gr$version-$packingfraction.tmp.sh";

    my $outfilename = "$dir/gr$version-$packingfraction.dat";
    my $command = "./radial-distribution-monte-carlo $packingfraction $acc $dr $outfilename";

    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
#SBATCH --output gr$version-$packingfraction.out

set -ev

hostname
date

echo I think this will take $memuse megs of memory
time nice -19 $command

cat > $outfilename.sum << ENDF
$sha1code
$command
ENDF

date

";
    close(SCRIPT);
    my $oldsha1 = "\n";
    if (-e "$outfilename.sum") {
      $oldsha1 = read_file("$outfilename.sum");
    }
    if ($oldsha1 eq "$sha1code\n$command\n") {
      print "Radial distribution with packing fraction $packingfraction is already done\n";
    } else {
      print "Submitting radial-distribution-monte-carlo with packing fraction $packingfraction\n";

      system("sbatch", $scriptname);
    }
}
