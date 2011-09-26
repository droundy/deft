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
	figs/n0.svg\
	figs/contact.svg\
	figs/free-energy.gp\
	figs/gHS-vs-n.gp\
	figs/sphere-radial-06-04.gp\
	figs/sphere-radial-08-13.gp\
	figs/box-100c--05,05,05-13.gp\
	figs/box-100s--05,05,05-13.gp\
	figs/box-110c--05,05,05-13.gp\
	figs/contact-08-013.pdf figs/contact-12-112.pdf figs/contact-16-265.pdf\

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

figs/gHS-vs-n.eps: figs/gHS-vs-n.dat

figs/free-energy.eps: figs/free-energy.dat

figs/sphere-radial-06-04.eps: figs/sphere-radial-06.0-04.dat
figs/sphere-radial-08-13.eps: figs/sphere-radial-08.0-13.dat figs/sphere-radial-12.0-112.dat \
	figs/sphere-radial-16.0-265.dat

figs/box-%.eps: figs/box-%.dat

figs/contact-08-013.pdf: figs/mc-04-013.dat figs/sphere-radial-08.0-13.dat figs/contact-density.py
	python figs/contact-density.py figs/mc-04-013.dat figs/sphere-radial-08.0-13.dat $@

figs/contact-12-112.pdf: figs/mc-06-112.dat figs/sphere-radial-12.0-112.dat figs/contact-density.py
	python figs/contact-density.py figs/mc-06-112.dat figs/sphere-radial-12.0-112.dat $@

figs/contact-16-265.pdf: figs/mc-08-265.dat figs/sphere-radial-16.0-265.dat figs/contact-density.py
	python figs/contact-density.py figs/mc-08-265.dat figs/sphere-radial-16.0-265.dat $@
