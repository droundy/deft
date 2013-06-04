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
all :  paper.pdf march-2013.pdf

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf

ALL_FIGURES=\
	figs/w2convolves.pdf figs/p-vs-packing.pdf figs/harmonic.pdf figs/w_3.pdf figs/w_2.pdf figs/w_1.pdf figs/step.pdf\
	figs/radial-distribution-10.pdf figs/walls-10.pdf

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

march-2013.pdf: march-2013.tex $(PDF_FIGURES)
	pdflatex $< pdflatex $< pdflatex $<

figs/w2convolves.pdf: figs/functions_plot.py
	python $<

figs/radial-distribution-10.pdf: figs/radial-distribution.py $(wildcard figs/mc-*.gradial)
	python $<

figs/walls-10.pdf: figs/plot-walls.py $(wildcard figs/mcwalls-*.dat)
	python $<

figs/p-vs-packing.pdf: figs/homogeneous.py $(wildcard figs/mc-*.prs)
	python $<

figs/harmonic.pdf: figs/harmonic.py
	python $^ $@

figs/w_3.pdf: figs/weight-functions.py
	python $^ $@

figs/w_2.pdf: figs/weight-functions.py
	python $^ $@

figs/w_1.pdf: figs/weight-functions.py
	python $^ $@

figs/step.pdf: figs/weight-functions.py
	python $^ $@
