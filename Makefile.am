.PHONY: paper project poster papers comparison slowwaterfigs haskell pdf

EXTRA_DIST = autogen.sh

NEW_CODE = src/new/test.cpp src/new/Minimize.cpp

FUNCTIONAL_CODE = src/lattice.cpp src/utilities.cpp \
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp \
	src/IdealGas.cpp src/ChemicalPotential.cpp \
	src/HardSpheres.cpp src/ExternalPotential.cpp \
	src/Functional.cpp \
	src/Gaussian.cpp src/Pow.cpp

GENERIC_CODE = $(FUNCTIONAL_CODE) \
	src/EffectivePotentialToDensity.cpp \
	src/equation-of-state.cpp src/water-constants.cpp \
	src/compute-surface-tension.cpp \
	src/WaterSaftFast.cpp src/WaterSaft_by_handFast.cpp \
	src/Minimizer.cpp src/Downhill.cpp \
	src/Precision.cpp src/ConjugateGradient.cpp \
	src/QuadraticLineMinimizer.cpp src/SteepestDescent.cpp

CONTACT_CODE = $(GENERIC_CODE) \
	src/ContactDensity.cpp

# src/HardSphereGasFast.cpp src/HardSphereGasRFFast.cpp \
#	src/HardSpheresFast.cpp src/HardSpheresWBFast.cpp src/HardSpheresWBm2Fast.cpp \
#	src/HardSpheresNoTensorFast.cpp

HASKELL_CODE = src/HardSpheresNoTensor2Fast.cpp \
	src/TensorWhiteBearFast.cpp src/WhiteBearMarkIIFast.cpp \
	src/TensorDensityXXFast.cpp src/VectorDensityXFast.cpp src/n2DensityFast.cpp \
	src/YuWuCorrelationFast.cpp src/SaftFluid2Fast.cpp \
	src/EntropySaftFluid2Fast.cpp src/CorrelationGrossCorrectFast.cpp \
	src/gSigmaSm2Fast.cpp src/gSigmaAm2Fast.cpp \
	src/gSigmaSFast.cpp src/gSigmaAFast.cpp \
	src/gSigmaA_automagicFast.cpp src/gSigmaA_by_handFast.cpp \
	src/SoftFluidFast.cpp src/HardFluidFast.cpp \
	src/WaterXFast.cpp src/HughesXFast.cpp \
	src/HughesHBFast.cpp

NEWHASKELL_CODE = src/new/WhiteBearFast.cpp
src/new/WhiteBearFast.h: src/new/WhiteBearFast.cpp
tests/new-hard-spheres.o: src/new/WhiteBearFast.h # doing this by hand is ugly!

ALL_CODE = $(GENERIC_CODE) $(HASKELL_CODE)

ALL_CONTACT_CODE = $(CONTACT_CODE) $(HASKELL_CODE)

bin_PROGRAMS = deft monte-carlo soft-monte-carlo pair-monte-carlo triplet-monte-carlo polyhedra-monte-carlo radial-distribution-monte-carlo # electrostatic-monte-carlo
deft_SOURCES = src/deft.cpp $(ALL_CODE)
tests_newcode_test_SOURCES = $(NEW_CODE)
monte_carlo_SOURCES = src/Monte-Carlo/monte-carlo.cpp \
	src/utilities.cpp
soft_monte_carlo_SOURCES = src/Monte-Carlo/soft-monte-carlo.cpp \
	src/utilities.cpp
pair_monte_carlo_SOURCES = src/Monte-Carlo/pair-monte-carlo.cpp \
	src/utilities.cpp
triplet_monte_carlo_SOURCES = src/Monte-Carlo/triplet-monte-carlo.cpp \
	src/utilities.cpp
polyhedra_monte_carlo_SOURCES = src/Monte-Carlo/polyhedra-monte-carlo.cpp \
	src/utilities.cpp
radial_distribution_monte_carlo_SOURCES = src/Monte-Carlo/radial-distribution-monte-carlo.cpp \
	src/utilities.cpp
electrostatic_monte_carlo_SOURCES = src/Monte-Carlo/electrostatic-monte-carlo.cpp \
	src/utilities.cpp
noinst_PROGRAMS =  \
	papers/water-saft/figs/surface-tension.mkdat \
	papers/water-saft/figs/equation-of-state.mkdat \
	papers/water-saft/figs/rods-in-water.mkdat \
	papers/water-saft/figs/four-rods-in-water.mkdat \
	papers/water-saft/figs/single-rod.mkdat \
	papers/water-saft/figs/hughes-single-rod.mkdat \
	papers/water-saft/figs/sphere.mkdat \
	papers/water-saft/figs/hughes-sphere.mkdat \
	papers/water-saft/figs/pressure-with-isotherms.mkdat \
	papers/hughes-saft/figs/surface-tension.mkdat \
	papers/hughes-saft/figs/equation-of-state.mkdat \
	papers/hughes-saft/figs/rods-in-water.mkdat \
	papers/hughes-saft/figs/four-rods-in-water.mkdat \
	papers/hughes-saft/figs/single-rod.mkdat \
	papers/hughes-saft/figs/single-rod-in-water-high-res.mkdat \
	papers/hughes-saft/figs/single-rod-in-water-low-res.mkdat \
	papers/hughes-saft/figs/sphere.mkdat \
	papers/hughes-saft/figs/pressure-with-isotherms.mkdat \
	papers/contact/figs/gHS-vs-n.mkdat \
	papers/contact/figs/free-energy.mkdat \
	papers/contact/figs/sphere.mkdat \
	papers/contact/figs/inner-sphere.mkdat \
	papers/contact/figs/test-particle-wall.mkdat \
	papers/contact/figs/walls.mkdat \
	papers/pair-correlation/figs/walls.mkdat \
	papers/pair-correlation/figs/sphere-with-wall.mkdat \
	papers/pair-correlation/figs/triplet-dft.mkdat \
	papers/fuzzy-fmt/figs/walls.mkdat

haskell: src/haskell/functionals.exe src/haskell/test.exe

# The following uses a weird trick.  We are using "pattern rule", to
# tell make that a single rule (calling another make, in this
# instance) can be used to generate several executables.  The pattern,
# in this case, matches the letter "s".
src/haskell/latex-functionals.%xe src/haskell/functionals.%xe src/haskell/newfunctionals.%xe src/haskell/test.%xe: $(wildcard src/haskell/*hs)
	cd src/haskell && $(MAKE) depend
	cd src/haskell && $(MAKE)

$(HASKELL_CODE): src/haskell/functionals.exe
	time src/haskell/functionals.exe $@ && echo timed $@

$(NEWHASKELL_CODE): src/haskell/newfunctionals.exe
	time src/haskell/newfunctionals.exe $@ && echo timed $@

tests/generated-haskell/nice-sum.h tests/generated-haskell/math.tex: src/haskell/test.exe
	test -d tests/generated-haskell || mkdir tests/generated-haskell
	rm -f tests/generated-haskell/*
	src/haskell/test.exe codegen
tests/generated-code.o: tests/generated-haskell/nice-sum.h

tests/new-generated-haskell/volume_minus_one_sqr.h: src/haskell/newfunctionals.exe
	test -d tests/new-generated-haskell || mkdir tests/new-generated-haskell
	rm -f tests/new-generated-haskell/*
	src/haskell/newfunctionals.exe tests
tests/new-generated-code.o: tests/new-generated-haskell/volume_minus_one_sqr.h

check_PROGRAMS = \
	tests/new-hard-spheres.test \
	tests/newcode.test \
	tests/generated-code.test \
	tests/new-generated-code.test \
	tests/saft.test \
	tests/memory.test \
	tests/eos.test \
	tests/print-iter.test \
	tests/ideal-gas.test \
	tests/convolve.test \
	tests/compare-gsigmas.test \
	tests/convolve-finite-difference.test \
	tests/precision.test \
	tests/functional-arithmetic.test \
	tests/functional-of-double.test \
	tests/surface-tension.test \
	tests/fftinverse.test \
	tests/eps.test

TESTS = tests/null.sh tests/print-iter.sh src/haskell/test.exe \
			$(check_PROGRAMS)
tests_saft_test_SOURCES = tests/saft.cpp $(ALL_CODE)
tests_generated_code_test_SOURCES = tests/generated-code.cpp $(GENERIC_CODE) \
	src/TensorWhiteBearFast.cpp src/gSigmaA_by_handFast.cpp src/gSigmaA_automagicFast.cpp
tests_new_generated_code_test_SOURCES = tests/new-generated-code.cpp $(GENERIC_CODE)
tests_new_hard_spheres_test_SOURCES = tests/new-hard-spheres.cpp $(GENERIC_CODE) $(NEWHASKELL_CODE)
tests_print_iter_test_SOURCES = tests/print-iter.cpp $(ALL_CODE)
tests_memory_test_SOURCES = tests/memory.cpp $(ALL_CONTACT_CODE)
tests_compare_gsigmas_test_SOURCES = tests/compare-gsigmas.cpp $(GENERIC_CODE) \
			src/gSigmaA_by_handFast.cpp src/gSigmaA_automagicFast.cpp
tests_convolve_test_SOURCES = tests/convolve.cpp $(GENERIC_CODE) \
			src/TensorDensityXXFast.cpp src/VectorDensityXFast.cpp src/n2DensityFast.cpp
tests_convolve_finite_difference_test_SOURCES = tests/convolve-finite-difference.cpp $(ALL_CODE)
tests_precision_test_SOURCES = tests/precision.cpp $(ALL_CODE)
tests_functional_arithmetic_test_SOURCES = tests/functional-arithmetic.cpp $(GENERIC_CODE)
tests_functional_of_double_test_SOURCES = tests/functional-of-double.cpp $(ALL_CONTACT_CODE)
tests_surface_tension_test_SOURCES = tests/surface-tension.cpp $(GENERIC_CODE)
tests_fftinverse_test_SOURCES = tests/fftinverse.cpp $(ALL_CODE)
tests_eps_test_SOURCES = tests/eps.cpp $(ALL_CODE)
tests_eos_test_SOURCES = tests/eos.cpp $(ALL_CODE)
tests_ideal_gas_test_SOURCES = tests/ideal-gas.cpp $(ALL_CODE)

papers_water_saft_figs_surface_tension_mkdat_SOURCES = papers/water-saft/figs/surface-tension.cpp $(ALL_CODE)
papers_water_saft_figs_equation_of_state_mkdat_SOURCES = papers/water-saft/figs/equation-of-state.cpp $(ALL_CODE)
papers_water_saft_figs_rods_in_water_mkdat_SOURCES = papers/water-saft/figs/rods-in-water.cpp $(ALL_CODE)
papers_water_saft_figs_four_rods_in_water_mkdat_SOURCES = papers/water-saft/figs/four-rods-in-water.cpp $(ALL_CODE)
papers_water_saft_figs_single_rod_mkdat_SOURCES = papers/water-saft/figs/single-rod.cpp $(ALL_CODE)
papers_water_saft_figs_hughes_single_rod_mkdat_SOURCES = papers/water-saft/figs/hughes-single-rod.cpp $(ALL_CODE)
papers_water_saft_figs_sphere_mkdat_SOURCES = papers/water-saft/figs/sphere.cpp $(ALL_CODE)
papers_water_saft_figs_hughes_sphere_mkdat_SOURCES = papers/water-saft/figs/hughes-sphere.cpp $(ALL_CODE)
papers_water_saft_figs_pressure_with_isotherms_mkdat_SOURCES = papers/water-saft/figs/pressure-with-isotherms.cpp $(ALL_CODE)

papers_hughes_saft_figs_surface_tension_mkdat_SOURCES = papers/hughes-saft/figs/surface-tension.cpp $(ALL_CODE)
papers_hughes_saft_figs_equation_of_state_mkdat_SOURCES = papers/hughes-saft/figs/equation-of-state.cpp $(ALL_CODE)
papers_hughes_saft_figs_rods_in_water_mkdat_SOURCES = papers/hughes-saft/figs/rods-in-water.cpp $(ALL_CODE)
papers_hughes_saft_figs_four_rods_in_water_mkdat_SOURCES = papers/hughes-saft/figs/four-rods-in-water.cpp $(ALL_CODE)
papers_hughes_saft_figs_single_rod_in_water_high_res_mkdat_SOURCES = papers/hughes-saft/figs/single-rod-in-water-high-res.cpp $(ALL_CODE)
papers_hughes_saft_figs_single_rod_in_water_low_res_mkdat_SOURCES = papers/hughes-saft/figs/single-rod-in-water-low-res.cpp $(ALL_CODE)
papers_hughes_saft_figs_single_rod_mkdat_SOURCES = papers/hughes-saft/figs/single-rod.cpp $(ALL_CODE)
papers_hughes_saft_figs_sphere_mkdat_SOURCES = papers/hughes-saft/figs/sphere.cpp $(ALL_CODE)
papers_hughes_saft_figs_pressure_with_isotherms_mkdat_SOURCES = papers/hughes-saft/figs/pressure-with-isotherms.cpp $(ALL_CODE)

slowwaterfigs: papers/hughes-saft/figs/single-rod-in-water-high-res.dat \
	$(wildcard papers/hughes-saft/figs/rods-in-water-*nm.dat) \
	$(wildcard papers/hughes-saft/figs/single-rod-*nm-energy.dat)

papers/hughes-saft/figs/rods-in-water-0%nm.dat: papers/hughes-saft/figs/rods-in-water.mkdat
	time $< $*

papers/hughes-saft/figs/single-rod-%nm-energy.dat: papers/hughes-saft/figs/single-rod.mkdat
	time $< $*

papers/water-saft/figs/rods-in-water-0%nm.dat: papers/water-saft/figs/rods-in-water.mkdat
	time $< $*

papers/water-saft/figs/single-rod-%nm-energy.dat: papers/water-saft/figs/single-rod.mkdat
	time $< $*

papers_contact_figs_gHS_vs_n_mkdat_SOURCES = papers/contact/figs/gHS-vs-n.cpp $(ALL_CONTACT_CODE)
papers_contact_figs_free_energy_mkdat_SOURCES = papers/contact/figs/free-energy.cpp $(CONTACT_CODE)
papers_contact_figs_sphere_mkdat_SOURCES = papers/contact/figs/sphere.cpp $(ALL_CONTACT_CODE)
papers_contact_figs_test_particle_wall_mkdat_SOURCES = papers/contact/figs/test-particle-wall.cpp $(ALL_CONTACT_CODE)
papers_contact_figs_inner_sphere_mkdat_SOURCES = papers/contact/figs/inner-sphere.cpp $(ALL_CONTACT_CODE)
papers_contact_figs_walls_mkdat_SOURCES = papers/contact/figs/walls.cpp $(ALL_CONTACT_CODE)

papers_pair_correlation_figs_walls_mkdat_SOURCES = papers/pair-correlation/figs/walls.cpp $(ALL_CONTACT_CODE)
papers_pair_correlation_figs_sphere_with_wall_mkdat_SOURCES = papers/pair-correlation/figs/sphere-with-wall.cpp $(ALL_CONTACT_CODE)
papers_pair_correlation_figs_triplet_dft_mkdat_SOURCES = papers/pair-correlation/figs/triplet-dft.cpp $(ALL_CONTACT_CODE)

papers_fuzzy_fmt_figs_walls_mkdat_SOURCES = papers/fuzzy-fmt/figs/walls.cpp $(ALL_CODE)

#noinst_PROGRAMS = optest

pdf: doc/Association.pdf doc/WhiteBear.pdf doc/TensorWhiteBear.pdf \
	doc/WhiteBearMarkII.pdf doc/Dispersion.pdf doc/SaftFluid.pdf \
	doc/SimpDispersion.pdf doc/EntropySaftFluid.pdf \
	doc/CorrelationS.pdf doc/CorrelationA.pdf doc/CorrelationAm2.pdf \
	doc/GradDispersion.pdf doc/JoinedGradDispersion.pdf doc/SimpGradDispersion.pdf

html: pdf papers
	jekyll build

clean-local:
	rm -rf html
	cd papers && $(MAKE) clean
	cd src/haskell && $(MAKE) clean

# Below, we remove ~/.matplotlib/fontList.cache, which can prevent
# crashes in matplotlib.  Hokey, but effective.  For a bit more
# information, see:
#
# https://github.com/matplotlib/matplotlib/issues/198

papers: all papers/contact/figs/gHS-vs-n.dat \
		papers/contact/figs/sphere.mkdat \
		papers/contact/figs/inner-sphere.mkdat \
		papers/contact/figs/test-particle-wall.mkdat \
		papers/contact/figs/walls.dat \
		papers/contact/figs/free-energy.dat \
		papers/pair-correlation/figs/walls.dat \
		papers/water-saft/figs/equation-of-state.dat \
		papers/water-saft/figs/surface-tension.dat \
		papers/water-saft/figs/pressure-with-isotherms.dat \
		papers/hughes-saft/figs/equation-of-state.dat \
		papers/hughes-saft/figs/surface-tension.dat \
		papers/hughes-saft/figs/pressure-with-isotherms.dat \
		papers/hughes-saft/figs/single-rod-in-water-low-res.dat
	rm -f ~/.matplotlib/fontList.cache # see above for explanation
	cd papers && $(MAKE)

papers/%.dat: papers/%.mkdat
	$<
doc/%.pdf: src/haskell/latex-functionals.exe
	$< $@
