#!/usr/bin/perl


open(SUBMIT, ">submit/paper.tex") or die "couldn't make submission paper $!";

open(PAPER, "paper.tex") or die "oops $!";
my $aminfigure = 0;
while (<PAPER>) {
  s/twocolumn/preprint/;
  if ($aminfigure) {
    if (/end{figure}/ or /end{figure\*/) {
      $aminfigure = 0;
    }
  } else {
    if (/begin{figure}/ or /begin{figure\*/) {
      $aminfigure = 1;
    } elsif (/end{document}/) {
      break;
    } else {
      print SUBMIT $_;
    }
  }
}
close(PAPER);

open(PAPER, "paper.tex") or die "oops $!";
$aminfigure = 0;
while (<PAPER>) {
  s/twocolumn/preprint/;
  if ($aminfigure) {
    if (/end{figure}/ or /end{figure\*/) {
      $aminfigure = 0;
    }
    print SUBMIT $_;
  } else {
    if (/begin{figure}/ or /begin{figure\*/) {
      $aminfigure = 1;
      print SUBMIT $_;
    }
  }
}
close(PAPER);

print SUBMIT "\\end{document}\n";
