#!/usr/bin/perl -w

use strict;

print "#########################################################################
# This makefile is automatically generated using find-deps.pl.  Please  #
# do not edit it!  If you think you need to edit it, please either edit #
# Makefile itself, or find-deps.pl.                                     #
#########################################################################
";

my %mainfiles;
my $makefile = ".make.depend1";

open(MAKEFILE, $makefile) or die "oops $!";
while (<MAKEFILE>) {
  if (/(.*).o : (.*hs)/) {
    my $fn = $1;
    if (open(HSFILE, "$2")) {
      my $line;
      foreach $line (<HSFILE>) {
        if ($line =~/^main/) {
          $mainfiles{$fn} = 1;
        }
      }
      close(HSFILE);
    }
  }
}
close(MAKEFILE);

my $main;
foreach $main (sort keys %mainfiles) {
  my %deps = ($main => 1);
  my $foundsomething = 1;

  while ($foundsomething) {
    $foundsomething = 0;
    open(MAKEFILE, $makefile) or die "oops $!";
    while (<MAKEFILE>) {
      if (/(.*).o : (.*).hi/) {
        #print "$1 depends on $2\n";
        if ($deps{$1} and not $deps{$2}) {
          #print "am looking for this!\n";
          $deps{$2} = 1;
          $foundsomething = 1
        } else {
          #print "no '$1'\n";
        }
      }
    }
    close(MAKEFILE);
  }

  print "\n$main.exe: ";
  my $name;
  foreach $name (sort keys %deps) {
    if ($name eq "") {
    } else {
      print "$name.o ";
    }
  }

  print "\n";
	print "\tghc \$(GHCFLAGS) -o \$@ \$^\n\n";
}
