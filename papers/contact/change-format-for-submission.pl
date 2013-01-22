#!/usr/bin/perl


open(SUBMIT, ">submit/paper.tex") or die "couldn't make submission paper $!";

open(PAPER, "paper.tex") or die "oops $!";
while (<PAPER>) {
  #s/twocolumn/preprint/;
  s/figs\///;
  print SUBMIT $_;
}
close(PAPER);

print SUBMIT "\\end{document}\n";
close(SUBMIT);
