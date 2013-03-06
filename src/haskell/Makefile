.SUFFIXES : .o .hs .hi .lhs

all: depend test.exe functionals.exe latex-functionals.exe newfunctionals.exe

clean:
	rm -f *.hi *.o test functionals latex-functionals newfunctionals

GHCFLAGS=-Wall -fno-warn-unused-do-bind -fspec-constr-count=20 -O2 -package process -package containers -package HUnit -Werror # -prof -auto-all

.o.hi:
	@:
.lhs.o:
	ghc -c $< $(GHCFLAGS)
.hs.o:
	ghc -c $< $(GHCFLAGS)

depend: make.depend

make.depend: find-deps.pl $(wildcard *hs)
	ghc -dep-makefile .make.depend1 -M $(GHCFLAGS) $(wildcard *hs)
	perl find-deps.pl > .make.depend2
	cat .make.depend2 .make.depend1 > make.depend

include make.depend
