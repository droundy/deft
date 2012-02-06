#!/usr/bin/perl

use strict;

my $iters  = 10000000;
my $acc = 0.001;
my $dir = "papers/contact/figs";
system "./monte-carlo 13 $iters $acc $dir/mc-04-013.dat outerSphere 4";
system "./monte-carlo 112 $iters $acc $dir/mc-06-112.dat outerSphere 6";
system "./monte-carlo 265 $iters $acc $dir/mc-08-265.dat outerSphere 8";
