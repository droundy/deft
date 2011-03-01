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
all :  paper.pdf

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf

ALL_FIGURES=\
	figs/correlation.gp \
	figs/equation-of-state.gp\
	figs/saturated-liquid-density.gp\
	figs/finding-vapor-pressure.gp\
	figs/near-critical-point.gp\
	figs/temperature-versus-density.gp\
	figs/surface-tension.gp\
	figs/surface-298.gp\
	figs/entropy.gp\

LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES)))))
FIGURES=$(patsubst %.jpg,%.eps, \
		$(patsubst %.png,%.eps,$(LESS_CANONICAL_FIGURES)))
PDF_FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))

paper.pdf: paper.tex paper.bib $(PDF_FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $<

figs/correlation.eps: figs/correlation.gp figs/gOO.exp figs/gOH.exp figs/gHH.exp

figs/finding-vapor-pressure.eps: figs/finding-vapor-pressure.dat
figs/equation-of-state.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
figs/saturated-liquid-density.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
figs/temperature-versus-density.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
figs/entropy.eps: figs/entropy.dat
figs/surface-tension.eps: figs/surface-tension.gp \
	figs/surface-tension.dat
figs/surface-298.eps: figs/surface-298.gp figs/surface-298.dat
figs/finding-vapor-pressure.eps: figs/finding-vapor-pressure.dat
figs/near-critical-point.eps: figs/finding-vapor-pressure.dat

