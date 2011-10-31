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
all :  paper.pdf poster.pdf

clean:
	rm -f *.aux *.dvi *.ps *.log *.pdf *.2602gf feynfile*.mf */*.eps */*.pdf

ALL_FIGURES=\
	figs/correlation.gp \
	figs/equation-of-state.gp\
	figs/saturated-liquid-density.gp\
	figs/finding-vapor-pressure.gp\
	figs/near-critical-point.gp\
	figs/temperature-versus-density.gp\
	figs/surface-tension.gp\
	figs/surface-298.gp\
	figs/entropy.gp\
	figs/entropy-at-690K.gp\
	figs/entropy-per-molecule.gp\
	figs/pressure-with-isotherms.gp\
	figs/cavitysize-060.gp\
	figs/cavity-cartoon.gp\
	figs/cavitysize-vs-energy-1D.gp\
	figs/xassoc-1D.gp\
	figs/entropy-1D.gp\
	figs/energy-1D.gp\
	figs/density-1D.gp\
	figs/sphere-energy-vs-diameter.gp\
	figs/sphere-density.gp\
	figs/energy-vs-diameter.gp\
	figs/density-rods-in-water.gp\
	figs/density-single-rod.gp\
	figs/density-vs-radius.gp\
	figs/rods-energy-vs-distance.gp\
	figs/energy-rods.gp\
	figs/osu-logo.svg\
	figs/xassoc-rods-in-water.gp\
	figs/xassoc-single-rod.gp\
	figs/two-rods-diagram.svg\
	figs/rods-diagram.svg\


LESS_CANONICAL_FIGURES= $(patsubst %.svg,%.eps, \
						$(patsubst %.gp,%.eps, \
						$(patsubst %.xgr,%.eps, \
						$(patsubst %.pl,%.eps, $(ALL_FIGURES)))))
FIGURES=$(patsubst %.eps,%.pdf, \
			$(LESS_CANONICAL_FIGURES))

project.pdf: project.tex paper.bib $(FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $< && pdflatex $<

paper.pdf: paper.tex paper.bib $(FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $<

poster.pdf: poster.tex paper.bib $(FIGURES)
	pdflatex $< && pdflatex $< && bibtex $* && pdflatex $< && pdflatex $<

figs/correlation.eps: figs/correlation.gp figs/gOO.exp figs/gOH.exp figs/gHH.exp

figs/finding-vapor-pressure.eps: figs/finding-vapor-pressure.dat
figs/equation-of-state.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
figs/saturated-liquid-density.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
figs/temperature-versus-density.eps: figs/equation-of-state.dat \
	figs/experimental-equation-of-state.dat
figs/entropy.eps: figs/entropy.dat
figs/entropy-at-690K.eps: figs/entropy-at-690K.dat
figs/entropy-per-molecule.eps: figs/entropy-per-molecule.gp figs/entropy.dat \
	figs/entropy-at-690K.dat 
figs/pressure-with-isotherms.eps: figs/pressure-with-isotherms.dat \
	figs/equation-of-state.dat
figs/surface-tension.eps: figs/surface-tension.gp \
	figs/surface-tension.dat
figs/surface-298.eps: figs/surface-298.gp figs/surface-298.dat
figs/finding-vapor-pressure.eps: figs/finding-vapor-pressure.dat
figs/near-critical-point.eps: figs/finding-vapor-pressure.dat
figs/cavitysize-060.eps: figs/cavitysize-06.0.dat
figs/cavity-cartoon.eps: figs/cavitysize-06.0.dat
figs/cavitysize-vs-energy-1D.eps: figs/constrained-water-1D.dat
figs/xassoc-1D.eps: figs/cavitysize-06.0.dat 
figs/entropy-1D.eps: figs/cavitysize-06.0.dat 
figs/energy-1D.eps: figs/cavitysize-06.0.dat
figs/density-1D.eps: figs/cavitysize-06.0.dat
figs/density-rods-in-water.eps: figs/rods-picture-01.0-00.2.dat \
	figs/rods-picture-01.0-00.4.dat
figs/xassoc-single-rod.gp: figs/single-rod-slice-00.6.dat
figs/xassoc-rods-in-water.gp: figs/rods-1nm-01.0.dat
figs/density-single-rod.eps: figs/single-rod-01.0.dat
figs/density-vs-radius.eps: figs/single-rod-slice-01.0.dat \
	figs/single-rod-res0.05-slice-01.0.dat figs/single-rod-res0.5-slice-01.0.dat 
figs/rods-energy-vs-distance.eps: figs/rods-in-water-00.6nm.dat \
	figs/rods-in-water-01.0nm.dat figs/rods-in-water-01.4nm.dat \
	figs/rods-in-water-01.8nm.dat figs/rods-in-water-02.0nm.dat
figs/energy-rods.eps: figs/rods-in-water-01.0nm.dat
figs/energy-vs-diameter.eps: figs/single-rod-in-water.dat

figs/sphere-energy-vs-diameter.eps: figs/sphere.dat
figs/sphere.dat: $(wildcard figs/sphere-%nm-energy.dat)
	echo 0 0 > $@
	cat figs/sphere-*nm-energy.dat >> $@
figs/sphere-density.eps: figs/rods-picture-01.0-00.2.dat \
	figs/rods-picture-01.0-00.4.dat
