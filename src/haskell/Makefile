all: test functionals

clean:
	rm -f *.hi *.o test

test: test.lhs $(wildcard *.lhs)
	ghc -Werror -Wall -fno-warn-unused-do-bind --make -o $@ $<

functionals: functionals.hs $(wildcard *.lhs)
	ghc -O2 -Wall -fno-warn-unused-do-bind --make -o $@ $<
