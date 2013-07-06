.SUFFIXES : .o .hs .hi .lhs

all: depend $(HASDIR)test.exe $(HASDIR)functionals.exe $(HASDIR)latex-functionals.exe $(HASDIR)newfunctionals.exe

clean:
	rm -f *.hi *.o *.exe

GHCFLAGS=-Wall -fno-warn-unused-do-bind -fspec-constr-count=20 -O2 -package process -package containers -package HUnit -Werror # -prof -auto-all

.o.hi:
	@:
$(HASDIR)%.o: %.hs
	ghc -c $< -o $@ $(GHCFLAGS) -odir $(HASDIR). -hidir $(HASDIR). -i$(HASDIR).

depend: find-deps.pl $(wildcard *hs)
	ghc -dep-makefile .make.depend1 -M $(GHCFLAGS) $(wildcard *hs)
	perl find-deps.pl > $(HASDIR)make.depend

include $(HASDIR)make.depend
