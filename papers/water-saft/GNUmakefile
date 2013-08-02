.SUFFIXES: .tex .dvi .ps .bib .bbl .pdf .fig .eps .aux .jpg .png .svg \
		.gp .py .mf .2602gf .pl .xgr
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
.py.pdf:
	python $<
.xgr.eps:
	gracebat -printfile $@ -hdevice EPS $<
.mf.2602gf:
	mf $<
all :  paper.pdf

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf
	rm -f paper/figs/*.tmp.sh

ALL_FIGURES=\
	figs/pressure-with-isotherms-truncated.gp\
	figs/surface-tension.gp\
	figs/energy-vs-diameter.gp\
	figs/density-rods-in-water.gp\
	figs/rods-energy-vs-distance.gp\
	figs/four-rods-energy-vs-distance.gp\
	figs/density-four-rods.gp\
	figs/sphere-energy-vs-diameter.gp\
	figs/density-sphere.gp\
	figs/density-compare.py\
	figs/single-rod-X-plot.py\
	figs/spheres-X-plot.py\
	figs/single-rod-broken-HB.py\
	figs/sphere-broken-HB.py



LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.py,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES))))))
FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))

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
figs/density-compare.pdf: figs/single-rod-slice-0.10.dat \
		figs/single-rod-slice-0.30.dat figs/single-rod-slice-0.60.dat \
		figs/single-rod-slice-1.00.dat figs/single-rod-slice-1.60.dat \
		figs/single-rod-slice-2.00.dat \
		figs/hughes-single-rod-slice-0.10.dat \
		figs/hughes-single-rod-slice-0.30.dat figs/hughes-single-rod-slice-0.60.dat \
		figs/hughes-single-rod-slice-1.00.dat figs/hughes-single-rod-slice-1.60.dat \
		figs/hughes-single-rod-slice-2.00.dat
figs/single-rod-X-plot.pdf: figs/single-rod-slice-0.10.dat \
		figs/single-rod-slice-0.30.dat figs/single-rod-slice-0.60.dat \
		figs/single-rod-slice-1.00.dat figs/single-rod-slice-1.60.dat \
		figs/single-rod-slice-2.00.dat figs/single-rod-slice-0.40.dat\
		figs/hughes-single-rod-slice-0.10.dat \
		figs/hughes-single-rod-slice-0.30.dat figs/hughes-single-rod-slice-0.60.dat \
		figs/hughes-single-rod-slice-1.00.dat figs/hughes-single-rod-slice-1.60.dat \
		figs/hughes-single-rod-slice-2.00.dat	figs/hughes-single-rod-slice-0.40.dat
figs/spheres-X-plot.pdf: figs/sphere-0.20.dat figs/sphere-0.40.dat \
		figs/sphere-0.80.dat figs/sphere-1.00.dat figs/sphere-1.20.dat \
		figs/sphere-1.40.dat figs/sphere-1.60.dat figs/sphere-2.00.dat \
		figs/hughes-sphere-0.20.dat figs/hughes-sphere-0.40.dat \
		figs/hughes-sphere-0.80.dat figs/hughes-sphere-1.00.dat figs/hughes-sphere-1.20.dat \
		figs/hughes-sphere-1.40.dat figs/hughes-sphere-1.60.dat figs/hughes-sphere-2.00.dat
figs/energy-vs-diameter.eps: figs/single-rod-in-water.dat
figs/single-rod-broken-HB.pdf: figs/single-rod-in-water.dat \
figs/hughes-single-rod-in-water.dat
figs/sphere-broken-HB.pdf: figs/hughes-sphere.dat figs/sphere.dat
singleRodEnergies := $(wildcard figs/single-rod-*nm-energy.dat)
figs/single-rod-in-water.dat: $(singleRodEnergies)
	cat $^ > $@
figs/single-rods-calculated-density.dat: figs/density_calc.py figs/single-rod-in-water.dat
	python figs/density_calc.py
figs/hughes-single-rod-in-water.dat: $(wildcard figs/hughes-single-rod-*nm-energy.dat)
	cat $^ > $@

figs/sphere-energy-vs-diameter.eps: figs/sphere.dat
sphereEnergies := $(wildcard figs/sphere-*nm-energy.dat)
figs/sphere.dat: $(sphereEnergies)
	cat $(sphereEnergies) > $@
figs/density-sphere.eps: $(wildcard figs/sphere-*.dat) $(wildcard figs/grsp*.dat)
hughesSphereEnergies := $(wildcard figs/hughes-sphere-*nm-energy.dat)
figs/hughes-sphere.dat: $(hughesSphereEnergies)
	cat $(hughesSphereEnergies) > $@
