.PHONY: paper project poster papers comparison slowfigs

EXTRA_DIST = autogen.sh

MONTECARLO_CODE = src/Monte-Carlo/monte-carlo.cpp \
	src/utilities.cpp \
	src/Monte-Carlo/fileio.cpp \
	src/Monte-Carlo/simulation.cpp 

FUNCTIONAL_CODE = src/lattice.cpp src/utilities.cpp \
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp \
	src/IdealGas.cpp src/ChemicalPotential.cpp src/GaussianPolynomial.cpp \
	src/HardSpheres.cpp src/ExternalPotential.cpp \
	src/Functional.cpp src/SAFT.cpp src/CallMe.cpp \
	src/Gaussian.cpp src/Pow.cpp src/Expression.cpp src/CSE.cpp

GENERIC_CODE = $(FUNCTIONAL_CODE) \
	src/ContactDensity.cpp \
	src/EffectivePotentialToDensity.cpp \
	src/equation-of-state.cpp src/water-constants.cpp \
	src/compute-surface-tension.cpp \
	src/Minimizer.cpp src/Downhill.cpp \
	src/Precision.cpp src/ConjugateGradient.cpp \
	src/QuadraticLineMinimizer.cpp src/SteepestDescent.cpp

OPTIMIZED_CODE = src/HardSphereGasFast.cpp src/HardSphereGasRFFast.cpp \
	src/SaftFluidFast.cpp \
	src/AssociationFast.cpp src/DispersionFast.cpp \
	src/HardSpheresFast.cpp src/HardSpheresRFFast.cpp \
	src/HardSpheresWBFast.cpp \
	src/HardSpheresTarazonaFast.cpp src/HardSpheresNoTensorFast.cpp

ALL_CODE = $(OPTIMIZED_CODE) $(GENERIC_CODE)

bin_PROGRAMS = deft monte-carlo analyze-monte-carlo
deft_SOURCES = src/deft.cpp $(ALL_CODE) 
monte_carlo_SOURCES = $(MONTECARLO_CODE)
analyze_monte_carlo_SOURCES = src/Monte-Carlo/analyze-monte-carlo.cpp \
	src/utilities.cpp
noinst_PROGRAMS = optimize-functionals tests/generating-code.test \
	paper/figs/surface-tension.mkdat \
	paper/figs/equation-of-state.mkdat \
	paper/figs/finding-vapor-pressure.mkdat \
	paper/figs/entropy.mkdat \
	paper/figs/entropy-at-690K.mkdat \
	paper/figs/constrained-water-1D.mkdat \
	paper/figs/constrained-water-2D.mkdat \
	paper/figs/rods-in-water.mkdat \
	paper/figs/four-rods-in-water.mkdat \
	paper/figs/single-rod-in-water.mkdat \
	paper/figs/single-rod-in-water-high-res.mkdat \
	paper/figs/single-rod-in-water-low-res.mkdat \
	paper/figs/sphere.mkdat \
	paper/figs/pressure-with-isotherms.mkdat \
	paper/figs/enthalpy.mkdat \
	papers/contact/figs/gHS-vs-n.mkdat \
	papers/contact/figs/free-energy.mkdat \
	papers/contact/figs/sphere.mkdat \
	papers/contact/figs/box.mkdat \
	comparisons/eta_effective_deft.test \
	comparisons/gHSetaeff_deft.test \
	comparisons/assoc_deft.test \
	comparisons/delta_deft.test \
	comparisons/x_deft.test \
	comparisons/gsw_deft.test \
	comparisons/da1dlam_deft.test \
	comparisons/da1deta_deft.test \
	comparisons/a2_sw_deft.test \
	comparisons/a1_sw_deft.test
optimize_functionals_SOURCES = src/optimize-functionals.cpp $(FUNCTIONAL_CODE)

$(OPTIMIZED_CODE): optimize-functionals
	./optimize-functionals $@

# The $(COMPARISONDATA) below is a fake dependency in order to ensure
# that the comparisons are run prior to the tests.
tests/generated/sum.h: tests/generating-code.test $(COMPARISONDATA)
	test -d tests/generated || mkdir tests/generated
	rm -f tests/generated/*
	tests/generating-code.test
tests/generated-code.o: tests/generated/sum.h

check_PROGRAMS = \
	tests/generated-code.test \
	tests/saft.test \
	tests/saft-dispersion.test \
	tests/optimized-code.test \
	tests/constrained-water-fast-slow.test \
	tests/expression.test \
	tests/functional-expressions.test \
	tests/print-iter.test \
	tests/memory.test \
	tests/ideal-gas.test \
	tests/hard-spheres-cavity.test \
	tests/hard-spheres-hard-wall.test \
	tests/hard-spheres.test \
	tests/convolve.test \
	tests/convolve-finite-difference.test \
	tests/simple-liquid.test \
	tests/precision.test \
	tests/functional-arithmetic.test \
	tests/functional-of-double.test \
	tests/surface-tension.test \
	tests/fftinverse.test \
	tests/eos.test \
	tests/eps.test

TESTS = tests/null.sh tests/print-iter.sh tests/comparisons.sh $(check_PROGRAMS)
tests_saft_test_SOURCES = tests/saft.cpp $(GENERIC_CODE)
tests_saft_dispersion_test_SOURCES = tests/saft-dispersion.cpp $(GENERIC_CODE)
tests_generated_code_test_SOURCES = tests/generated-code.cpp $(GENERIC_CODE)
tests_optimized_code_test_SOURCES = tests/optimized-code.cpp $(ALL_CODE)
tests_constrained_water_fast_slow_test_SOURCES = tests/constrained-water-fast-slow.cpp $(ALL_CODE)
tests_generating_code_test_SOURCES = tests/generating-code.cpp $(GENERIC_CODE)
tests_functional_expressions_test_SOURCES = tests/functional-expressions.cpp $(GENERIC_CODE)
tests_expression_test_SOURCES = tests/expression.cpp $(GENERIC_CODE)
tests_print_iter_test_SOURCES = tests/print-iter.cpp $(ALL_CODE)
tests_memory_test_SOURCES = tests/memory.cpp $(ALL_CODE)
tests_hard_spheres_cavity_test_SOURCES = tests/hard-spheres-cavity.cpp $(ALL_CODE)
tests_hard_spheres_hard_wall_test_SOURCES = tests/hard-spheres-hard-wall.cpp $(ALL_CODE)
tests_hard_spheres_test_SOURCES = tests/hard-spheres.cpp $(ALL_CODE)
tests_convolve_test_SOURCES = tests/convolve.cpp $(ALL_CODE)
tests_convolve_finite_difference_test_SOURCES = tests/convolve-finite-difference.cpp $(ALL_CODE)
tests_precision_test_SOURCES = tests/precision.cpp $(ALL_CODE)
tests_functional_arithmetic_test_SOURCES = tests/functional-arithmetic.cpp $(GENERIC_CODE)
tests_functional_of_double_test_SOURCES = tests/functional-of-double.cpp $(GENERIC_CODE)
tests_surface_tension_test_SOURCES = tests/surface-tension.cpp $(ALL_CODE)
tests_simple_liquid_test_SOURCES = tests/simple-liquid.cpp $(ALL_CODE)
tests_fftinverse_test_SOURCES = tests/fftinverse.cpp $(ALL_CODE)
tests_eps_test_SOURCES = tests/eps.cpp $(ALL_CODE)
tests_eos_test_SOURCES = tests/eos.cpp $(ALL_CODE)
tests_ideal_gas_test_SOURCES = tests/ideal-gas.cpp $(ALL_CODE)

paper_figs_surface_tension_mkdat_SOURCES = paper/figs/surface-tension.cpp $(ALL_CODE)
paper_figs_equation_of_state_mkdat_SOURCES = paper/figs/equation-of-state.cpp $(ALL_CODE)
paper_figs_finding_vapor_pressure_mkdat_SOURCES = paper/figs/finding-vapor-pressure.cpp $(ALL_CODE)
paper_figs_entropy_mkdat_SOURCES = paper/figs/entropy.cpp $(GENERIC_CODE)
paper_figs_entropy_at_690K_mkdat_SOURCES = paper/figs/entropy-at-690K.cpp $(GENERIC_CODE)
paper_figs_constrained_water_1D_mkdat_SOURCES = paper/figs/constrained-water-1D.cpp $(ALL_CODE)
paper_figs_constrained_water_2D_mkdat_SOURCES = paper/figs/constrained-water-2D.cpp $(ALL_CODE)
paper_figs_rods_in_water_mkdat_SOURCES = paper/figs/rods-in-water.cpp $(ALL_CODE)
paper_figs_four_rods_in_water_mkdat_SOURCES = paper/figs/four-rods-in-water.cpp $(ALL_CODE)
paper_figs_single_rod_in_water_high_res_mkdat_SOURCES = paper/figs/single-rod-in-water-high-res.cpp $(ALL_CODE)
paper_figs_single_rod_in_water_low_res_mkdat_SOURCES = paper/figs/single-rod-in-water-low-res.cpp $(ALL_CODE)
paper_figs_single_rod_in_water_mkdat_SOURCES = paper/figs/single-rod-in-water.cpp $(ALL_CODE)
paper_figs_sphere_mkdat_SOURCES = paper/figs/sphere.cpp $(ALL_CODE)
paper_figs_pressure_with_isotherms_mkdat_SOURCES = paper/figs/pressure-with-isotherms.cpp $(ALL_CODE)
paper_figs_enthalpy_mkdat_SOURCES = paper/figs/enthalpy.cpp $(ALL_CODE)

papers_contact_figs_gHS_vs_n_mkdat_SOURCES = papers/contact/figs/gHS-vs-n.cpp $(GENERIC_CODE)
papers_contact_figs_free_energy_mkdat_SOURCES = papers/contact/figs/free-energy.cpp $(GENERIC_CODE)
papers_contact_figs_sphere_mkdat_SOURCES = papers/contact/figs/sphere.cpp $(ALL_CODE)
papers_contact_figs_box_mkdat_SOURCES = papers/contact/figs/box.cpp $(ALL_CODE)

comparisons_eta_effective_deft_test_SOURCES = comparisons/eta_effective_deft.cpp $(GENERIC_CODE)
comparisons_assoc_deft_test_SOURCES = comparisons/assoc_deft.cpp $(GENERIC_CODE)
comparisons_delta_deft_test_SOURCES = comparisons/delta_deft.cpp $(GENERIC_CODE)
comparisons_x_deft_test_SOURCES = comparisons/x_deft.cpp $(GENERIC_CODE)
comparisons_gsw_deft_test_SOURCES = comparisons/gsw_deft.cpp $(GENERIC_CODE)
comparisons_da1dlam_deft_test_SOURCES = comparisons/da1dlam_deft.cpp $(GENERIC_CODE)
comparisons_da1deta_deft_test_SOURCES = comparisons/da1deta_deft.cpp $(GENERIC_CODE)
comparisons_a1_sw_deft_test_SOURCES = comparisons/a1_sw_deft.cpp $(GENERIC_CODE)
comparisons_a2_sw_deft_test_SOURCES = comparisons/a2_sw_deft.cpp $(GENERIC_CODE)
comparisons_gHSetaeff_deft_test_SOURCES = comparisons/gHSetaeff_deft.cpp $(GENERIC_CODE)

#noinst_PROGRAMS = optest

#BUILT_SOURCES = patchlevel.h

#patchlevel.h: update_patchlevel.sh @DARCS_INVENTORY@ $(DFT_SOURCES) $(CMND_SOURCES)
#	sh update_patchlevel.sh
# Build static html docs suitable for being shipped in the software
# package. This depends on ikiwiki being installed to build the docs.

# ifeq ($(shell which ikiwiki),)
# IKIWIKI=echo "** ikiwiki not found" >&2 ; echo ikiwiki
# else
# IKIWIKI=ikiwiki
# endif

html:
	ikiwiki `pwd` html -v --wikiname Deft --plugin=goodstuff \
		--exclude=html --exclude=Makefile.am

clean-local:
	rm -rf .ikiwiki html
	cd paper && $(MAKE) clean
	cd papers && $(MAKE) clean
	rm -f comparisons/*.o comparisons/*.mkdat

paper: paper/figs/equation-of-state.dat \
	paper/figs/finding-vapor-pressure.dat \
	paper/figs/surface-tension.dat \
	paper/figs/enthalpy.dat \
	paper/figs/entropy.dat \
	paper/figs/entropy-at-690K.dat \
	paper/figs/constrained-water-1D.dat \
	paper/figs/constrained-water-2D.dat \
	paper/figs/pressure-with-isotherms.dat \
	paper/figs/single-rod-in-water-low-res.dat
	cd paper && $(MAKE)

papers: papers/contact/figs/gHS-vs-n.dat \
	papers/contact/figs/sphere.mkdat papers/contact/figs/box.mkdat \
	papers/contact/figs/free-energy.dat
	cd papers && $(MAKE)

project: paper
	cd paper && $(MAKE) project.pdf

poster: paper
	cd paper && $(MAKE) poster.pdf

slowfigs: paper/figs/rods-in-water.dat \
	paper/figs/four-rods-in-water.dat\
	paper/figs/single-rod-in-water.dat \
	paper/figs/single-rod-in-water-high-res.dat \
	paper/figs/sphere.dat \
	papers/contact/figs/mc-06-112.dat \
	papers/contact/figs/box.dat

paper/figs/%.dat: paper/figs/%.mkdat
	$<
papers/contact/figs/%.dat: papers/contact/figs/%.mkdat
	$<

COMPARISONDATA = comparisons/eta_effective_vrpack.dat comparisons/eta_effective_deft.dat \
	comparisons/gHSetaeff_vrpack.dat comparisons/gHSetaeff_deft.dat \
	comparisons/assoc_vrpack.dat comparisons/assoc_deft.dat \
	comparisons/delta_vrpack.dat comparisons/delta_deft.dat \
	comparisons/x_vrpack.dat comparisons/x_deft.dat \
	comparisons/gsw_vrpack.dat comparisons/gsw_deft.dat \
	comparisons/da1dlam_vrpack.dat comparisons/da1dlam_deft.dat \
	comparisons/da1deta_vrpack.dat comparisons/da1deta_deft.dat \
	comparisons/a1_sw_vrpack.dat comparisons/a1_sw_deft.dat \
	comparisons/a2_sw_vrpack.dat comparisons/a2_sw_deft.dat

comparisons/eta_effective_vrpack.test: comparisons/eta_effective_vrpack.f fortran/vrpack.f
	gfortran -o $@ $<
# The following rule works with any of the comparisons...
comparisons/%.dat: comparisons/%.test
	$< > $@

comparisons/a1_sw_vrpack.test: comparisons/a1_sw_vrpack.f fortran/vrpack.f
	gfortran -o $@ $<
comparisons/a2_sw_vrpack.test: comparisons/a2_sw_vrpack.f fortran/vrpack.f
	gfortran -o $@ $<
comparisons/gHSetaeff_vrpack.test: comparisons/gHSetaeff_vrpack.f fortran/vrpack.f
	gfortran -o $@ $<
comparisons/assoc_vrpack.test: comparisons/assoc_vrpack.f fortran/saft_assoc.f fortran/vrpack.f fortran/saft_assoc.f
	gfortran -o $@ $<
comparisons/delta_vrpack.test: comparisons/delta_vrpack.f fortran/saft_assoc.f fortran/vrpack.f fortran/saft_assoc.f
	gfortran -o $@ $<
comparisons/x_vrpack.test: comparisons/x_vrpack.f fortran/saft_assoc.f fortran/vrpack.f fortran/saft_assoc.f
	gfortran -o $@ $<
comparisons/gsw_vrpack.test: comparisons/gsw_vrpack.f fortran/vrpack.f
	gfortran -o $@ $<
comparisons/da1dlam_vrpack.test: comparisons/da1dlam_vrpack.f fortran/vrpack.f
	gfortran -o $@ $<
comparisons/da1deta_vrpack.test: comparisons/da1deta_vrpack.f fortran/vrpack.f
	gfortran -o $@ $<

comparison: $(COMPARISONDATA)
	tests/comparisons.sh

MCITERATIONS = 1000000
papers/contact/figs/mc-06-112.dat: monte-carlo
	time ./monte-carlo 6 112 $(MCITERATIONS) $@


papers/contact/figs/mc-density-%.dat: papers/contact/figs/mc-%.dat \
		analyze-monte-carlo
	time ./analyze-monte-carlo $@ $<
