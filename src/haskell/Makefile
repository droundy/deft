all: test

clean:
	rm -f *.hi *.o test

test: test.lhs $(wildcard *.lhs)
	ghc -Werror -Wall -fno-warn-unused-do-bind --make -o test test.lhs
