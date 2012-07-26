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
all :  paper.pdf march-2012.pdf

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf
	rm -f paper/figs/*.tmp.sh

ALL_FIGURES=\
	figs/correlation.gp \
	figs/equation-of-state.gp\
	figs/saturated-liquid-density.gp\
	figs/temperature-versus-density.gp\
	figs/surface-tension.gp\
	figs/surface-298.gp\
	figs/pressure-with-isotherms.gp\
	figs/pressure-with-isotherms-truncated.gp\
	figs/sphere-energy-vs-diameter.gp\
	figs/density-sphere.gp\
	figs/energy-vs-diameter.gp\
	figs/density-rods-in-water.gp\
	figs/density-four-rods.gp\
	figs/density-single-rod.gp\
	figs/rods-energy-vs-distance.gp\
	figs/four-rods-energy-vs-distance.gp\
	figs/density-four-rods.gp\
	figs/osu-logo.svg\
	figs/stadium-vs-circles.svg\
	figs/two-rods-diagram.svg\
	figs/four-rods-diagram.svg\
	figs/rods-diagram.svg\


LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES)))))
FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))

paper.pdf: paper.tex paper.bib $(FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $< && pdflatex $<

march-2012.pdf: march-2012.tex $(FIGURES)
	pdflatex $< && pdflatex $<

figs/correlation.eps: figs/correlation.gp figs/gOO.exp figs/gOH.exp figs/gHH.exp

figs/equation-of-state.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
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
figs/density-rods-in-water.eps: figs/rods-slice-01.0-00.5.dat figs/rods-slice-01.0-00.6.dat
figs/density-four-rods.eps: figs/four-rods-01.0-1.22.dat \
	figs/four-rods-01.0-1.37.dat
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
figs/density-sphere.eps: $(wildcard figs/sphere-*.dat)
