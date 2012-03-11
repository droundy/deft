all: test functionals

clean:
	rm -f *.hi *.o test

GHCFLAGS=-Werror -Wall -fno-warn-unused-do-bind -O2

test: test.lhs $(wildcard *hs)
	ghc $(GHCFLAGS) --make -o $@ $<

# I make functionals depend on test, so they won't get confused by
# compiling both files at once (which can mean compiling CodeGen.lhs
# twice at the same time).
functionals: functionals.hs test $(wildcard *hs)
	ghc $(GHCFLAGS) --make -o $@ $<
