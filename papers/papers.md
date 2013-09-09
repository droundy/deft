# Deft papers

Here you can access several papers that present results obtained
using the Deft code.  The code that generated the results for each of
these papers is available in the Deft repository, and can serve as an
example of how to use Deft.

"Improved association in a classical density functional theory for
  water", E. J. Krebs, J. B. Sculte and D. Roundy, *submitted to
  J. Chem. Phys.* [(pdf)](water-saft/paper.pdf)

"A classical density-functional theory for describing water
interfaces", J. Hughes, E. J. Krebs, D. Roundy, *J. Chem. Phys.*
(2013). [(pdf)](hughes-saft/paper.pdf)

"Using fundamental measure theory to treat the correlation function of
the inhomogeneous hard-sphere fluid", J. B. Schulte, P. A. Kreitzberg,
C. V. Haglund, D. Roundy, *Phys. Rev. E*
(2013). [(pdf)](contact/paper.pdf)

You may also be interested in the
[slides](../talks/colloquium/slides.pdf) for a colloquium introducing
the main results of the above three papers.


# Information for Deft developers

The `papers/` directory is where we have the latex code and code to
generate the figures that will produce papers that are describe deft
and the methods it implements.

Each paper is in a subdirectory, which contains all the code required
to create the paper.  These subdirectories are themselves listed in
`papers/GNUmakefile`.

Within each paper directory is a `GNUmakefile` that controls the
building of the paper itself.  There is a `figs/` subdirectory, which
contains all the code (and data) required to generate the figures for
that paper.

