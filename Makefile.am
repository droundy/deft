EXTRA_DIST = autogen.sh

GENERIC_CODE = src/lattice.cpp src/utilities.cpp \
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp \
	src/IdealGas.cpp src/ChemicalPotential.cpp \
	src/HardSpheres.cpp \
	src/Minimizer.cpp src/Downhill.cpp src/ExternalPotential.cpp \
	src/Precision.cpp src/ConjugateGradient.cpp \
	src/QuadraticLineMinimizer.cpp src/SteepestDescent.cpp src/GaussianPolynomial.cpp \
	src/Gaussian.cpp src/Pow.cpp src/Functional.cpp src/Expression.cpp \
	src/EffectivePotentialToDensity.cpp

bin_PROGRAMS = deft optimize-functionals
deft_SOURCES = src/deft.cpp $(GENERIC_CODE) src/HardSpheresFast.h
optimize_functionals_SOURCES = src/optimize-functionals.cpp $(GENERIC_CODE)

src/HardSpheresFast.h: optimize-functionals
	./optimize-functionals

tests/generated/sum.h: tests/generating-code.test
	test -d tests/generated || mkdir tests/generated
	tests/generating-code.test
tests/generated-code.o: tests/generated/sum.h
tests/hard-spheres-fast.o: src/HardSpheresFast.h

check_PROGRAMS = \
	tests/generating-code.test \
	tests/generated-code.test \
	tests/timings.test \
	tests/expression.test \
	tests/functional-expressions.test \
	tests/print-iter.test \
	tests/memory.test \
	tests/ideal-gas.test \
	tests/hard-spheres-fast.test \
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
tests_generating_code_test_SOURCES = tests/generating-code.cpp $(GENERIC_CODE)
tests_timings_test_SOURCES = tests/timings.cpp $(GENERIC_CODE)
tests_functional_expressions_test_SOURCES = tests/functional-expressions.cpp $(GENERIC_CODE)
tests_expression_test_SOURCES = tests/expression.cpp $(GENERIC_CODE)
tests_print_iter_test_SOURCES = tests/print-iter.cpp $(GENERIC_CODE)
tests_memory_test_SOURCES = tests/memory.cpp $(GENERIC_CODE)
tests_hard_spheres_fast_test_SOURCES = tests/hard-spheres-fast.cpp $(GENERIC_CODE)
tests_hard_spheres_cavity_test_SOURCES = tests/hard-spheres-cavity.cpp $(GENERIC_CODE)
tests_hard_spheres_hard_wall_test_SOURCES = tests/hard-spheres-hard-wall.cpp $(GENERIC_CODE)
tests_hard_spheres_test_SOURCES = tests/hard-spheres.cpp $(GENERIC_CODE)
tests_convolve_test_SOURCES = tests/convolve.cpp $(GENERIC_CODE)
tests_convolve_finite_difference_test_SOURCES = tests/convolve-finite-difference.cpp $(GENERIC_CODE)
tests_precision_test_SOURCES = tests/precision.cpp $(GENERIC_CODE)
tests_functional_arithmetic_test_SOURCES = tests/functional-arithmetic.cpp $(GENERIC_CODE)
tests_functional_of_double_test_SOURCES = tests/functional-of-double.cpp $(GENERIC_CODE)
tests_surface_tension_test_SOURCES = tests/surface-tension.cpp $(GENERIC_CODE)
tests_simple_liquid_test_SOURCES = tests/simple-liquid.cpp $(GENERIC_CODE)
tests_fftinverse_test_SOURCES = tests/fftinverse.cpp $(GENERIC_CODE)
tests_eps_test_SOURCES = tests/eps.cpp $(GENERIC_CODE)
tests_ideal_gas_test_SOURCES = tests/ideal-gas.cpp $(GENERIC_CODE)
tests_ideal_gas_at_pressure_test_SOURCES = tests/ideal-gas-at-pressure.cpp $(GENERIC_CODE)
tests_ideal_gas_optimization_test_SOURCES = tests/ideal-gas-optimization.cpp $(GENERIC_CODE)
tests_ideal_gas_external_potential_test_SOURCES = tests/ideal-gas-external-potential.cpp $(GENERIC_CODE)

#noinst_PROGRAMS = optest
#optest_SOURCES = optest.cpp $(GENERIC_CODE)

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
