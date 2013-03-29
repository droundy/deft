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

ALL_FIGURES=\

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

arxiv/figs/%: figs/%
	cp $< $@

arxiv.tar.gz: paper.pdf jeff-notes.pdf$(ARXIV_FIGURES)
	cp paper.tex paper.bbl jeff-notes.tex jeff-notes.bbl arxiv/
	tar cvzf arxiv.tar.gz arxiv
