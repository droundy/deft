.SUFFIXES : .o .hs .hi .lhs

all: depend test functionals latex-functionals newfunctionals

clean:
	rm -f *.hi *.o test functionals latex-functionals newfunctionals

GHCFLAGS=-Wall -fno-warn-unused-do-bind -fspec-constr-count=20 -O2 -package process -package containers -package HUnit # -prof -auto-all

.o.hi:
	@:
.lhs.o:
	ghc -c $< $(GHCFLAGS)
.hs.o:
	ghc -c $< $(GHCFLAGS)

CODEGEN_O=Expression.o CodeGen.o Statement.o LatexDouble.o
NEWCODEGEN_O=Expression.o NewCode.o Statement.o LatexDouble.o

test: test.o SomeFunctionals.o WhiteBear.o FMT.o $(CODEGEN_O)
	$(MAKE) depend
	ghc $(GHCFLAGS) -o $@ $^

functionals: functionals.o SomeFunctionals.o WhiteBear.o FMT.o Quantum.o $(CODEGEN_O)
	ghc $(GHCFLAGS) -o $@ $^

newfunctionals: newfunctionals.o WhiteBear.o FMT.o $(NEWCODEGEN_O)
	ghc $(GHCFLAGS) -o $@ $^

latex-functionals: latex-functionals.o SomeFunctionals.o WhiteBear.o FMT.o Latex.o Quantum.o $(CODEGEN_O)
	ghc $(GHCFLAGS) -o $@ $^

depend:
	ghc -M $(GHCFLAGS) $(wildcard *hs)

# DO NOT DELETE: Beginning of Haskell dependencies
LatexDouble.o : LatexDouble.hs
Expression.o : Expression.lhs
Expression.o : LatexDouble.hi
FMT.o : FMT.hs
FMT.o : Expression.hi
WhiteBear.o : WhiteBear.hs
WhiteBear.o : FMT.hi
WhiteBear.o : Expression.hi
SomeFunctionals.o : SomeFunctionals.hs
SomeFunctionals.o : Expression.hi
SomeFunctionals.o : WhiteBear.hi
SomeFunctionals.o : FMT.hi
Statement.o : Statement.lhs
Statement.o : Expression.hi
Latex.o : Latex.hs
Latex.o : Statement.hi
Latex.o : Expression.hi
NewCode.o : NewCode.lhs
NewCode.o : Expression.hi
NewCode.o : Statement.hi
newfunctionals.o : newfunctionals.hs
newfunctionals.o : WhiteBear.hi
newfunctionals.o : NewCode.hi
CodeGen.o : CodeGen.lhs
CodeGen.o : Expression.hi
CodeGen.o : Statement.hi
Quantum.o : Quantum.hs
Quantum.o : CodeGen.hi
latex-functionals.o : latex-functionals.hs
latex-functionals.o : Expression.hi
latex-functionals.o : Latex.hi
latex-functionals.o : Quantum.hi
latex-functionals.o : WhiteBear.hi
latex-functionals.o : SomeFunctionals.hi
functionals.o : functionals.hs
functionals.o : Quantum.hi
functionals.o : WhiteBear.hi
functionals.o : SomeFunctionals.hi
functionals.o : CodeGen.hi
test.o : test.lhs
test.o : FMT.hi
test.o : WhiteBear.hi
test.o : SomeFunctionals.hi
test.o : CodeGen.hi
# DO NOT DELETE: End of Haskell dependencies
