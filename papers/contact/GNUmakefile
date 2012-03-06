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
	figs/free-energy.gp\
	figs/gHS-vs-n.gp\
	figs/contact-08-013.pdf figs/contact-12-112.pdf figs/contact-16-265.pdf\
	figs/walls-10.pdf figs/walls-30.pdf figs/walls-40.pdf \
	figs/inner-4-30.pdf figs/inner-4-40.pdf figs/inner-8-30.pdf \

LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES)))))
FIGURES=$(patsubst %.jpg,%.eps, \
		$(patsubst %.png,%.eps,$(LESS_CANONICAL_FIGURES)))
PDF_FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))

paper.pdf: paper.tex formulas.tex paper.bib $(PDF_FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $< && pdflatex $<

figs/gHS-vs-n.eps: figs/gHS-vs-n.dat

figs/free-energy.eps: figs/free-energy.dat

#figs/box-%.eps: figs/box-%.dat

figs/contact-08-013.pdf: figs/contact-density.py figs/mc-04-013.dat figs/sphereWB-radial-08.0-13.dat figs/sphereWBT-radial-08.0-13.dat figs/sphereWBm2-radial-08.0-13.dat
	python $^ $@

figs/contact-12-112.pdf: figs/contact-density.py figs/mc-06-112.dat figs/sphereWB-radial-12.0-112.dat figs/sphereWBT-radial-12.0-112.dat figs/sphereWBm2-radial-12.0-112.dat
	python $^ $@

figs/contact-16-265.pdf: figs/contact-density.py figs/mc-08-265.dat figs/sphereWB-radial-16.0-265.dat figs/sphereWBT-radial-16.0-265.dat figs/sphereWBm2-radial-16.0-265.dat
	python $^ $@

figs/walls-10.pdf: figs/plot-walls.py figs/mc-walls-15-300.dat figs/wallsWB-0.10.dat figs/wallsWBT-0.10.dat figs/wallsWBm2-0.10.dat
	python $^ $@

figs/walls-30.pdf: figs/plot-walls.py figs/mc-walls-15-300.dat figs/wallsWB-0.30.dat figs/wallsWBT-0.30.dat figs/wallsWBm2-0.30.dat
	python $^ $@

figs/walls-40.pdf: figs/plot-walls.py figs/mc-walls-15-300.dat figs/wallsWB-0.40.dat figs/wallsWBT-0.40.dat figs/wallsWBm2-0.40.dat
	python $^ $@

figs/inner-4-30.pdf: figs/plot-inner.py figs/mc-innerSphere-20-02-760.dat figs/inner-sphereWB-04.0-0.30.dat figs/inner-sphereWBT-04.0-0.30.dat figs/inner-sphereWBm2-04.0-0.30.dat
	python $^ $@

figs/inner-4-40.pdf: figs/plot-inner.py figs/mc-innerSphere-20-02-760.dat figs/inner-sphereWB-04.0-0.40.dat figs/inner-sphereWBT-04.0-0.40.dat figs/inner-sphereWBm2-04.0-0.40.dat
	python $^ $@

figs/inner-8-30.pdf: figs/plot-inner.py figs/mc-innerSphere-20-02-760.dat figs/inner-sphereWB-08.0-0.30.dat figs/inner-sphereWBT-08.0-0.30.dat figs/inner-sphereWBm2-08.0-0.30.dat
	python $^ $@
