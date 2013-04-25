#!/usr/bin/perl -w

use strict;
use Digest::SHA1 qw(sha1_hex);
use File::Slurp qw(read_file);

my $dir = "papers/fuzzy-fmt/figs";
my $sha1code = sha1_hex(read_file("src/Monte-Carlo/soft-monte-carlo.cpp"));

system "make soft-monte-carlo";

my $iters  = 99999999999999;
my $acc = 0.001;


my $dd;
foreach $dd ([20,382,1],[20,573,1],[20,764,1],[20,954,1]) {
    my $radius = sprintf("%02.0f", $$dd[0]);
    my $N = sprintf("%03.0f", $$dd[1]);
    my $kT = sprintf("%03.5f", $$dd[2]);

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 0.001*($N) + 30; # It's a very hokey guess

    my $scriptname = "papers/fuzzy-fmt/figs/mc-soft-sphere-homogeneous-$radius-$N-$kT.tmp.sh";

    my $outfilename = "$dir/mc-soft-homogenous-$radius-$N-$kT.dat";
    my $command = "./soft-monte-carlo $N $iters $acc $outfilename periodxyz $radius kT $kT";

    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user jeffschulte83\@gmail.com
#SBATCH --output mc-soft-hogeneous-$radius-$N.out

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
      print "Monte-carlo in cubic cavity with inner sphere radius  has already been done.\n"
    } else {
      print "Submitting monte-carlo in cubic cavity with inner sphere radius \n";

      system("sbatch", $scriptname);
    }
}

foreach $dd ([20,382,1],[20,573,1],[20,764,1],[20,954,1]){
    my $len = sprintf("%02.0f", $$dd[0]);
    my $N = sprintf("%03.0f", $$dd[1]);
    my $kT = sprintf("%03.5f", $$dd[2]);

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 0.001*($N) + 30; # It's a very hokey guess

    my $scriptname = "papers/fuzzy-fmt/figs/mc-soft-walls-$len-$N-$kT.tmp.sh";

    my $outfilename = "$dir/mc-soft-walls-$len-$N-$kT.dat";
    my $command = "./soft-monte-carlo $N $iters $acc $outfilename periodx $len periody $len wallz $len flatdiv kT $kT";

    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user jeffschulte83\@gmail.com
#SBATCH --output mc-soft-walls-$len-$N-$kT.out

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
      print "Monte-carlo in cubic cavity with length $len wallz kT $kT and flat divisions has already been done.\n"
    } else {
      print "Submitting monte-carlo in cubic cavity with length $len wallz kT $kT and flat divisions\n";

      system("sbatch", $scriptname);
    }
}

