.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .mf .2602gf .pl .xgr
.tex.pdf:
	pdflatex $<
.ps.pdf:
	ps2pdf $<
.jpg.eps:
	convert $< $@
.pl.eps:
	perl $<
.png.eps:
	convert $< $@
.eps.pdf:
	epstopdf $<
.svg.eps:
	inkscape --export-eps $@ $<
.dvi.ps:
	dvips -t letter -Ppdf -o $@ $<
.bib.bbl:
	bibtex $*
.fig.eps:
	fig2dev -L eps $< $@

.gp.eps:
	gnuplot $<
.xgr.eps:
	gracebat -printfile $@ -hdevice EPS $<
.mf.2602gf:
	mf $<
all :  project.pdf presentation.pdf

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf
	rm -f *.bbl *.blg *.toc *.vrb *.snm *.nav

ALL_FIGURES=\
	figs/equation-of-state.gp\
	figs/saturated-liquid-density.gp\
	figs/finding-vapor-pressure.gp\
	figs/near-critical-point.gp\
	figs/temperature-versus-density.gp\
	figs/surface-tension.gp\
	figs/surface-298.gp\
	figs/entropy.gp\
	figs/pressure-with-isotherms.gp\
	figs/xassoc-1D.gp\
	figs/energy-1D.gp\
	figs/density-1D.gp\
	figs/energy-vs-diameter.gp\
	figs/density-rods-in-water.gp\
	figs/density-single-rod.gp\
	figs/density-vs-radius.gp\
	figs/rods-energy-vs-distance.gp\
	figs/energy-rods.gp\
	figs/osu-logo.svg\
	figs/xassoc-single-rod.gp\
	figs/two-rods-diagram.svg\
	figs/single-rod-diagram.svg\
	figs/rods-diagram.svg\
	figs/hydrogen-bonding.svg\
	figs/slit-diagram.svg\


LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES)))))
FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))

project.pdf: project.tex paper.bib $(FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $< && pdflatex $<

presentation.pdf: presentation.tex $(FIGURES)
	pdflatex $< && pdflatex $<
