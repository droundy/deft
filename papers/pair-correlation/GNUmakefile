.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .mf .2602gf .pl .xgr
.tex.pdf:
	pdflatex $<
#.ps.pdf:
#	ps2pdf $<
.jpg.eps:
	convert $< $@
.pl.eps:
	perl $<
.png.eps:
	convert $< $@
.eps.pdf:
	epstopdf $<
.pdf.ps:
	pdftops $< $@
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
all :  paper.pdf jeff-notes.pdf

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf
	rm -f arxiv/*
	mkdir arxiv/figs

ALL_FIGURES= figs/pair-correlation-2-395-0.pdf figs/pair-correlation-2-395-5.pdf \
	figs/pair-correlation-2-395-10.pdf \
	figs/dadz-3-2.pdf figs/dadz-3-3.pdf figs/dadz-2-2.pdf figs/dadz-2-3.pdf\
	figs/pair-correlation-pretty-1.pdf \
	figs/pair-correlation-pretty-2.pdf \
	figs/pair-correlation-pretty-3.pdf \
	figs/pair-correlation-pretty-4.pdf \
	figs/ghs-g.pdf \

LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES)))))
FIGURES=$(patsubst %.pdf,%.ps, \
		$(patsubst %.jpg,%.eps, \
		$(patsubst %.png,%.eps,$(LESS_CANONICAL_FIGURES))))
ARXIV_FIGURES=$(patsubst figs/%,arxiv/figs/%,$(FIGURES))
PDF_FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))

paper.pdf: paper.tex paper.bib jeff-notes.tex jeff-notes.bib $(PDF_FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $< && pdflatex $<

figs/ghs-g.pdf: figs/plot-ghs.py
	python figs/plot-ghs.py

figs/dadz-3-2.pdf: figs/plot-da_dz.py
	python figs/plot-da_dz.py .3 2.01

figs/dadz-3-3.pdf: figs/plot-da_dz.py
	python figs/plot-da_dz.py .3 3.0

figs/dadz-2-2.pdf: figs/plot-da_dz.py
	python figs/plot-da_dz.py .2 2.01

figs/dadz-2-3.pdf: figs/plot-da_dz.py
	python figs/plot-da_dz.py .2 3.0

figs/pair-correlation-2-395-0.pdf: figs/plot-mc-modified.py
	python figs/plot-mc-modified.py .2 3.95 0.0

figs/pair-correlation-2-395-5.pdf: figs/plot-mc-modified.py
	python figs/plot-mc-modified.py .2 3.95 0.50

figs/pair-correlation-2-395-10.pdf: figs/plot-mc-modified.py
	python figs/plot-mc-modified.py .2 3.95 1.0

figs/pair-correlation-pretty-1.pdf: figs/plot-path-and-2d.py $(wildcard figs/*.dat)
	python figs/plot-path-and-2d.py 0.1
figs/pair-correlation-pretty-2.pdf: figs/plot-path-and-2d.py $(wildcard figs/*.dat)
	python figs/plot-path-and-2d.py 0.2
figs/pair-correlation-pretty-3.pdf: figs/plot-path-and-2d.py $(wildcard figs/*.dat)
	python figs/plot-path-and-2d.py 0.3
figs/pair-correlation-pretty-4.pdf: figs/plot-path-and-2d.py $(wildcard figs/*.dat)
	python figs/plot-path-and-2d.py 0.4



arxiv/figs/%: figs/%
	cp $< $@

arxiv.tar.gz: paper.pdf jeff-notes.pdf$(ARXIV_FIGURES)
	cp paper.tex paper.bbl jeff-notes.tex jeff-notes.bbl arxiv/
	tar cvzf arxiv.tar.gz arxiv
