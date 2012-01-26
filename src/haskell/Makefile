all: test functionals

clean:
	rm -f *.hi *.o test

GHCFLAGS=-Werror -Wall -fno-warn-unused-do-bind

test: test.lhs $(wildcard *.lhs)
	ghc $(GHCFLAGS) --make -o $@ $<

functionals: functionals.hs $(wildcard *.lhs)
	ghc $(GHCFLAGS) --make -o $@ $<
