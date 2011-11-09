.PHONY: test

all: test

clean:
	rm -f *.hi *.o

test: test.lhs
	ghc --make -o test test.lhs
