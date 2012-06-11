.SUFFIXES : .o .hs .hi .lhs

all: depend test functionals latex-functionals

clean:
	rm -f *.hi *.o test functionals latex-functionals

GHCFLAGS=-Werror -Wall -fno-warn-unused-do-bind -O2 -package containers -package HUnit # -prof -auto-all
GHCFLAGS=-Wall -fno-warn-unused-do-bind -O2 -package containers -package HUnit # -prof -auto-all

.o.hi:
	@:
.lhs.o:
	ghc -c $< $(GHCFLAGS)
.hs.o:
	ghc -c $< $(GHCFLAGS)

CODEGEN_O=Expression.o CodeGen.o Statement.o LatexDouble.o

test: test.o SomeFunctionals.o $(CODEGEN_O)
	$(MAKE) depend
	ghc $(GHCFLAGS) -o $@ $^

functionals: functionals.o SomeFunctionals.o Quantum.o $(CODEGEN_O)
	ghc $(GHCFLAGS) -o $@ $^

latex-functionals: latex-functionals.o SomeFunctionals.o Latex.o Quantum.o $(CODEGEN_O)
	ghc $(GHCFLAGS) -o $@ $^

depend:
	ghc -M $(GHCFLAGS) $(wildcard *hs)

# DO NOT DELETE: Beginning of Haskell dependencies
LatexDouble.o : LatexDouble.hs
Expression.o : Expression.lhs
Expression.o : LatexDouble.hi
Statement.o : Statement.lhs
Statement.o : Expression.hi
Latex.o : Latex.hs
Latex.o : Statement.hi
Latex.o : Expression.hi
CodeGen.o : CodeGen.lhs
CodeGen.o : Expression.hi
CodeGen.o : Statement.hi
SomeFunctionals.o : SomeFunctionals.hs
SomeFunctionals.o : CodeGen.hi
Quantum.o : Quantum.hs
Quantum.o : CodeGen.hi
latex-functionals.o : latex-functionals.hs
latex-functionals.o : Expression.hi
latex-functionals.o : Latex.hi
latex-functionals.o : Quantum.hi
latex-functionals.o : SomeFunctionals.hi
functionals.o : functionals.hs
functionals.o : Quantum.hi
functionals.o : SomeFunctionals.hi
functionals.o : CodeGen.hi
test.o : test.lhs
test.o : SomeFunctionals.hi
test.o : CodeGen.hi
# DO NOT DELETE: End of Haskell dependencies
