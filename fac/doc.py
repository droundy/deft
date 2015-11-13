#!/usr/bin/python3

import facfile

doc = facfile.facfile('.doc.fac')

for pdf in """ Association WhiteBear TensorWhiteBear WhiteBearMarkII Dispersion SaftFluid
               SimpDispersion EntropySaftFluid GradDispersion JoinedGradDispersion
               SW_liquid
               SimpGradDispersion SFMT """.split():
    doc.default('src/haskell/latex-functionals.exe doc/%s.pdf' % pdf,
                ['src/haskell/latex-functionals.exe'],
                ['doc/%s.pdf' % pdf])
