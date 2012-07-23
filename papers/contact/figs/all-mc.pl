#!/usr/bin/perl -w

use strict;
use Digest::SHA1 qw(sha1_hex);
use File::Slurp qw(read_file);

my $dir = "papers/contact/figs";
my $sha1code = sha1_hex(read_file("src/Monte-Carlo/monte-carlo.cpp"));

system "make monte-carlo";

my $iters  = 3*1000*1000*1000*1000;
my $acc = 0.001;

my $dd;
foreach $dd ([20,832], [20,1718], [20,2629], [20,3543], [20,4475],
    [16,430], [16,891], [16,1370], [16,1880], [16,2375],
    [12,185], [12,386], [12,598], [12,822], [12,995]) {
    my $radius = sprintf("%02.0f", $$dd[0]);
    my $N = sprintf("%03.0f", $$dd[1]);

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 0.001*($N) + 30; # It's a very hokey guess

    my $scriptname = "papers/contact/figs/mc-sphere-$radius-$N.tmp.sh";

    my $outfilename = "$dir/mc-$radius-$N.dat";
    my $command = "./monte-carlo $N $iters $acc $outfilename outerSphere $radius";

    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user jeffschulte83\@gmail.com
#SBATCH --output mc-sphere-$radius-$N.out

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
      print "Monte-carlo in cavity with radius $radius has already been done.\n"
    } else {
      print "Submitting monte-carlo in cavity with radius $radius\n";

      system("sbatch", $scriptname);
    }
}

foreach $dd ([20,195], [20,390], [20,605], [20,816], [20,1020],[24,1320],[24,1650]) {
    my $len = sprintf("%02.0f", $$dd[0]);
    my $N = sprintf("%03.0f", $$dd[1]);

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 0.001*($N) + 30; # It's a very hokey guess

    my $scriptname = "papers/contact/figs/mc-walls-$len-$N.tmp.sh";

    my $outfilename = "$dir/mc-walls-$len-$N.dat";
    my $command = "./monte-carlo $N $iters $acc $outfilename periodx $len periody $len wallz $len flatdiv";

    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user jeffschulte83\@gmail.com
#SBATCH --output mc-walls-$len-$N.out

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
      print "Monte-carlo in cubic cavity with length $len wallz and flat divisions has already been done.\n"
    } else {
      print "Submitting monte-carlo in cubic cavity with length $len wallz and flat divisions\n";

      system("sbatch", $scriptname);
    }
}

#[2,24,1492], [2,25,1119], [4,28,1601]){

foreach $dd ([2,24,501], [2,24,751], [2,24,1153], [6,32,1400]){ #, [2,24,1488], [2,25,1120], [4,28,1602], [6,32,1602]){   #[2,24,330],[2,24,494], [2,25,1118], [2,24,1491], [2,24,659], [2,24,824], [2,24,990], [2,24,1320], [2,24,1485], [4,28,550], [4,28,750], [4,28,850], [4,28,950], [4,28,1100], [4,28,1300],[4,28,1400],[4,28,1500],[4,28,1600], [6,32,1400], [6,32,1500],[6,32,1600]){

    my $innerRad = sprintf("%02.0f", $$dd[0]);
    my $len = sprintf("%02.0f", $$dd[1]);
    my $N = sprintf("%03.0f", $$dd[2]);

    # Here I estimate the amount of memory that will be needed...
    my $memuse = sprintf "%.0f", 0.001*($N) + 30; # It's a very hokey guess

    my $scriptname = "papers/contact/figs/mc-periodic-innerSphere-$len-$innerRad-$N.tmp.sh";

    my $outfilename = "$dir/mc-innerSphere-$len-$innerRad-$N.dat";
    my $command = "./monte-carlo $N $iters $acc $outfilename innerSphere $innerRad periodx $len periody $len periodz $len";

    open SCRIPT, ">$scriptname" or die $!;
    print SCRIPT "#!/bin/sh
#SBATCH --mem-per-cpu=$memuse
##SBATCH --mail-type ALL
##SBATCH --mail-user jeffschulte83\@gmail.com
#SBATCH --output mc-innerSphere-$len-$innerRad-$N.out

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
      print "Monte-carlo in cubic cavity with inner sphere radius $innerRad has already been done.\n"
    } else {
      print "Submitting monte-carlo in cubic cavity with inner sphere radius $innerRad\n";

      system("sbatch", $scriptname);
    }
}

