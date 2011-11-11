.PHONY: test

all: test

clean:
	rm -f *.hi *.o

test: test.lhs
	ghc -Werror -Wall --make -o test test.lhs
