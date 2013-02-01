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
	figs/density-08-013.pdf figs/density-12-112.pdf figs/density-16-265.pdf figs/w2convolves.pdf\

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

figs/w2convolves.pdf: figs/functions_plot.py
	python $^ $@

figs/density-08-013.pdf: figs/density.py figs/mc-04-013.dat figs/sphereWB-radial-08.0-13.dat figs/sphereWBT-radial-08.0-13.dat
	python $^ $@

figs/density-12-112.pdf: figs/density.py figs/mc-06-112.dat figs/sphereWB-radial-12.0-112.dat figs/sphereWBT-radial-12.0-112.dat
	python $^ $@

figs/density-16-265.pdf: figs/density.py figs/mc-08-265.dat figs/sphereWB-radial-16.0-265.dat figs/sphereWBT-radial-16.0-265.dat
	python $^ $@
