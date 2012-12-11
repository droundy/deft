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
all :  paper.pdf arxiv.tar.gz submit.tar.gz

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf
	rm -f arxiv/*
	mkdir arxiv/figs
	rm -rf submit submit.pdf

ALL_FIGURES=\
	figs/n0.svg figs/nA.svg \
	figs/contact.svg\
	figs/free-energy.gp\
	figs/gHS-vs-n.gp\
	figs/gnn-sphere-08-013.pdf figs/gnn-sphere-12-112.pdf figs/gnn-sphere-16-265.pdf\
	figs/outer-16-04.pdf \
	figs/outer-16-01.pdf \
	figs/walls-10.pdf figs/walls-30.pdf figs/walls-40.pdf \
	figs/inner-4-10.pdf figs/inner-4-40.pdf \
	figs/inner-12-10.pdf figs/inner-12-40.pdf \
	figs/test-particle-10.pdf figs/test-particle-40.pdf \

LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES)))))
FIGURES=$(patsubst %.pdf,%.ps, \
		$(patsubst %.jpg,%.eps, \
		$(patsubst %.png,%.eps,$(LESS_CANONICAL_FIGURES))))
PDF_FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))
ARXIV_FIGURES=$(patsubst figs/%,arxiv/figs/%,$(PDF_FIGURES))
SUBMIT_FIGURES=$(patsubst figs/%,submit/figs/%, $(PDF_FIGURES))

paper.pdf: paper.tex paper.bib $(PDF_FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $< && pdflatex $<

figs/gHS-vs-n.eps: figs/gHS-vs-n.dat

figs/free-energy.eps: figs/free-energy.dat

#figs/box-%.eps: figs/box-%.dat

figs/outer-20-0%.pdf: figs/plot-outer.py figs/outer-sphereWB-20-00.%.dat
	python $< 20 $* $@
figs/outer-16-0%.pdf: figs/plot-outer.py figs/outer-sphereWB-16-00.%.dat
	python $< 16 $* $@
figs/outer-12-0%.pdf: figs/plot-outer.py figs/outer-sphereWB-12-00.%.dat
	python $< 12 $* $@

figs/gnn-sphere-08-013.pdf: figs/plot-gnn-sphere.py figs/mc-04-013.dat figs/sphereWB-radial-08.0-13.dat figs/sphereWBT-radial-08.0-13.dat figs/sphereWBm2-radial-08.0-13.dat
	python $^ $@

figs/gnn-sphere-12-112.pdf: figs/plot-gnn-sphere.py figs/mc-06-112.dat figs/sphereWB-radial-12.0-112.dat figs/sphereWBT-radial-12.0-112.dat figs/sphereWBm2-radial-12.0-112.dat
	python $^ $@

figs/gnn-sphere-16-265.pdf: figs/plot-gnn-sphere.py figs/mc-08-265.dat figs/sphereWB-radial-16.0-265.dat figs/sphereWBT-radial-16.0-265.dat figs/sphereWBm2-radial-16.0-265.dat
	python $^ $@

figs/walls-10.pdf: figs/plot-walls.py figs/mc-walls-20-196.dat figs/wallsWB-0.10.dat figs/wallsWBT-0.10.dat figs/wallsWBm2-0.10.dat
	python $^ $@

figs/walls-30.pdf: figs/plot-walls.py figs/mc-walls-20-605.dat figs/wallsWB-0.30.dat figs/wallsWBT-0.30.dat figs/wallsWBm2-0.30.dat
	python $^ $@

figs/walls-40.pdf: figs/plot-walls.py figs/mc-walls-20-817.dat figs/wallsWB-0.40.dat figs/wallsWBT-0.40.dat figs/wallsWBm2-0.40.dat
	python $^ $@

figs/gnn-walls-10.pdf: figs/plot-gnn-walls.py figs/mc-walls-20-193.dat figs/wallsWB-0.10.dat figs/wallsWBT-0.10.dat figs/wallsWBm2-0.10.dat
	python $^ $@

figs/gnn-walls-30.pdf: figs/plot-gnn-walls.py figs/mc-walls-20-605.dat figs/wallsWB-0.30.dat figs/wallsWBT-0.30.dat figs/wallsWBm2-0.30.dat
	python $^ $@

figs/gnn-walls-40.pdf: figs/plot-gnn-walls.py figs/mc-walls-20-816.dat figs/wallsWB-0.40.dat figs/wallsWBT-0.40.dat figs/wallsWBm2-0.40.dat
	python $^ $@

figs/inner-4-%0.pdf: figs/plot-inner.py figs/inner-sphereWB-04.0-0.%0.dat $(wildcard figs/mc-inner-sphere*.dat)
	python $< 2 $* $@

figs/inner-12-%0.pdf: figs/plot-inner.py figs/inner-sphereWB-12.0-0.%0.dat $(wildcard figs/mc-inner-sphere*.dat)
	python $< 6 $* $@

figs/ghs.pdf: figs/plot-ghs.py $(wildcard figs/mc-inner-sphere*.dat)
	python $< $@

test-particle-0.10 := $(wildcard figs/test-particle-cont-dens-*-0.10.dat)
figs/test-particle-cont-dens-0.10.dat: $(test-particle-0.10)
	cat $^ > $@

test-particle-0.40 := $(wildcard figs/test-particle-cont-dens-*-0.40.dat)
figs/test-particle-cont-dens-0.40.dat: $(test-particle-0.40)
	cat $^ > $@

figs/test-particle-10.pdf: figs/plot-test-particle-cont-dens.py figs/mc-walls-20-195.dat figs/wallsWB-0.10.dat figs/wallsWBT-0.10.dat figs/wallsWBm2-0.10.dat figs/test-particle-cont-dens-0.10.dat
	python $^ $@

figs/test-particle-40.pdf: figs/plot-test-particle-cont-dens.py figs/mc-walls-20-816.dat figs/wallsWB-0.40.dat figs/wallsWBT-0.40.dat figs/wallsWBm2-0.40.dat figs/test-particle-cont-dens-0.40.dat
	python $^ $@

arxiv/figs/%: figs/%
	mkdir -p arxiv/figs
	cp $< $@

arxiv.tar.gz: paper.pdf $(ARXIV_FIGURES)
	cp paper.tex paper.bbl arxiv/
	tar cvzf arxiv.tar.gz arxiv

submit/figs/%: figs/%
	mkdir -p submit/figs
	cp $< $@

submit.tar.gz: paper.pdf change-format-for-submission.pl $(SUBMIT_FIGURES)
	perl change-format-for-submission.pl
	cp paper.bbl submit/
	tar cvzf submit.tar.gz submit
	cd submit && pdflatex paper && pdflatex paper
	mv submit/paper.pdf submit.pdf
