EXTRA_DIST = autogen.sh

FUNCTIONAL_CODE = src/lattice.cpp src/utilities.cpp \
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp \
	src/IdealGas.cpp src/ChemicalPotential.cpp src/GaussianPolynomial.cpp \
	src/HardSpheres.cpp src/ExternalPotential.cpp src/Functional.cpp \
	src/SAFT.cpp \
	src/EffectivePotentialToDensity.cpp \
	src/Gaussian.cpp src/Pow.cpp src/Expression.cpp src/CSE.cpp

GENERIC_CODE = $(FUNCTIONAL_CODE) \
	src/Minimizer.cpp src/Downhill.cpp \
	src/Precision.cpp src/ConjugateGradient.cpp \
	src/QuadraticLineMinimizer.cpp src/SteepestDescent.cpp

OPTIMIZED_CODE = src/HardSphereGasFast.cpp src/HardSphereGasRFFast.cpp \
	src/SaftFluidFast.cpp \
	src/HardSpheresFast.cpp src/HardSpheresRFFast.cpp \
	src/HardSpheresTarazonaFast.cpp src/HardSpheresNoTensor.cpp

ALL_CODE = $(OPTIMIZED_CODE) $(GENERIC_CODE)

bin_PROGRAMS = deft
deft_SOURCES = src/deft.cpp $(ALL_CODE)
noinst_PROGRAMS = optimize-functionals
optimize_functionals_SOURCES = src/optimize-functionals.cpp $(FUNCTIONAL_CODE)

$(OPTIMIZED_CODE): optimize-functionals
	./optimize-functionals $@

tests/generated/sum.h: tests/generating-code.test
	test -d tests/generated || mkdir tests/generated
	tests/generating-code.test
tests/generated-code.o: tests/generated/sum.h

check_PROGRAMS = \
	tests/generating-code.test \
	tests/generated-code.test \
	tests/optimized-code.test \
	tests/timings.test \
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
	tests/ideal-gas-at-pressure.test \
	tests/ideal-gas-optimization.test \
	tests/ideal-gas-external-potential.test \
	tests/eps.test

TESTS = tests/null.sh tests/print-iter.sh $(check_PROGRAMS)
tests_generated_code_test_SOURCES = tests/generated-code.cpp $(GENERIC_CODE)
tests_optimized_code_test_SOURCES = tests/optimized-code.cpp $(ALL_CODE)
tests_generating_code_test_SOURCES = tests/generating-code.cpp $(GENERIC_CODE)
tests_timings_test_SOURCES = tests/timings.cpp $(ALL_CODE)
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
tests_functional_of_double_test_SOURCES = tests/functional-of-double.cpp $(ALL_CODE)
tests_surface_tension_test_SOURCES = tests/surface-tension.cpp $(ALL_CODE)
tests_simple_liquid_test_SOURCES = tests/simple-liquid.cpp $(ALL_CODE)
tests_fftinverse_test_SOURCES = tests/fftinverse.cpp $(ALL_CODE)
tests_eps_test_SOURCES = tests/eps.cpp $(ALL_CODE)
tests_ideal_gas_test_SOURCES = tests/ideal-gas.cpp $(ALL_CODE)
tests_ideal_gas_at_pressure_test_SOURCES = tests/ideal-gas-at-pressure.cpp $(ALL_CODE)
tests_ideal_gas_optimization_test_SOURCES = tests/ideal-gas-optimization.cpp $(ALL_CODE)
tests_ideal_gas_external_potential_test_SOURCES = tests/ideal-gas-external-potential.cpp $(ALL_CODE)

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

local-clean:
	rm -rf .ikiwiki html
