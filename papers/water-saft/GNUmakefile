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
all :  paper.pdf arxiv.tar.gz submit.tar.gz

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf
	rm -f paper/figs/*.tmp.sh
	rm -rf arxiv arxiv.tar.gz
	rm -rf submit submit.tar.gz

ALL_FIGURES=\
	figs/pressure-with-isotherms-truncated.gp\
	figs/surface-tension.gp\
	figs/energy-vs-diameter.gp\
	figs/density-single-rod.gp\
	figs/density-rods-in-water.gp\
	figs/rods-energy-vs-distance.gp\
	figs/four-rods-energy-vs-distance.gp\
	figs/density-four-rods.gp\
	figs/sphere-energy-vs-diameter.gp\
	figs/density-sphere.gp\
	figs/single-rod-X-plot.gp



LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES)))))
FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))
ARXIV_FIGURES=$(patsubst figs/%,arxiv/figs/%, \
	$(patsubst %.pdf,%.ps, $(LESS_CANONICAL_FIGURES)))
SUBMIT_FIGURES=$(patsubst figs/%,submit/figs/%, $(LESS_CANONICAL_FIGURES))

paper.pdf: paper.tex paper.bib $(FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $< && pdflatex $<

figs/correlation.eps: figs/correlation.gp figs/gOO.exp figs/gOH.exp figs/gHH.exp

figs/saturated-liquid-density.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
figs/temperature-versus-density.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
figs/pressure-with-isotherms.eps: figs/pressure-with-isotherms.dat \
	figs/equation-of-state.dat
figs/pressure-with-isotherms-truncated.eps: figs/pressure-with-isotherms.dat \
	figs/equation-of-state.dat
figs/surface-tension.eps: figs/surface-tension.gp \
	figs/surface-tension.dat
figs/surface-298.eps: figs/surface-298.gp figs/surface-298.dat
figs/density-rods-in-water.eps: figs/rods-slice-01.2-00.6.dat figs/rods-slice-01.2-00.7.dat
figs/density-four-rods.eps: figs/four-rods-01.2-1.44.dat \
	figs/four-rods-01.2-1.47.dat
figs/rods-energy-vs-distance.eps: figs/rods-in-water-00.4nm.dat \
	figs/rods-in-water-00.8nm.dat figs/rods-in-water-01.2nm.dat \
	figs/rods-in-water-01.6nm.dat figs/rods-in-water-02.0nm.dat \
	figs/rods-in-water-02.4nm.dat
figs/four-rods-energy-vs-distance.eps: figs/four-rods-in-water-00.4nm.dat \
	figs/four-rods-in-water-00.8nm.dat figs/four-rods-in-water-01.2nm.dat \
	figs/four-rods-in-water-01.6nm.dat figs/four-rods-in-water-02.0nm.dat \
	figs/four-rods-in-water-02.4nm.dat
figs/density-single-rod.eps: figs/single-rod-slice-00.1.dat \
		figs/single-rod-slice-00.3.dat figs/single-rod-slice-00.6.dat \
		figs/single-rod-slice-01.0.dat figs/single-rod-slice-01.6.dat \
		figs/single-rod-slice-02.0.dat figs/single-rod-slice-00.6.dat \
		figs/single-rods-calculated-density.dat
figs/single-rod-X-plot.eps: figs/single-rod-slice-00.1.dat \
		figs/single-rod-slice-00.3.dat figs/single-rod-slice-00.6.dat \
		figs/single-rod-slice-01.0.dat figs/single-rod-slice-01.6.dat \
		figs/single-rod-slice-02.0.dat figs/single-rod-slice-00.6.dat \
		figs/hughes-single-rod-slice-00.1.dat \
		figs/hughes-single-rod-slice-00.3.dat figs/hughes-single-rod-slice-00.6.dat \
		figs/hughes-single-rod-slice-01.0.dat figs/hughes-single-rod-slice-01.6.dat \
		figs/hughes-single-rod-slice-02.0.dat figs/hughes-single-rod-slice-00.6.dat
figs/energy-vs-diameter.eps: figs/single-rod-in-water.dat
singleRodEnergies := $(wildcard figs/single-rod-*nm-energy.dat)
figs/single-rod-in-water.dat: $(singleRodEnergies)
	cat $^ > $@
figs/single-rods-calculated-density.dat: figs/density_calc.py figs/single-rod-in-water.dat
	python figs/density_calc.py

figs/sphere-energy-vs-diameter.eps: figs/sphere.dat
sphereEnergies := $(wildcard figs/sphere-*nm-energy.dat)
figs/sphere.dat: $(sphereEnergies)
	cat $(sphereEnergies) > $@
figs/density-sphere.eps: $(wildcard figs/sphere-*.dat) $(wildcard figs/grsp*.dat)

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
	zip -r submit submit
	cd submit && latex paper && latex paper
	mv submit/paper.dvi submit.dvi
