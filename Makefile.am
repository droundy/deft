EXTRA_DIST = autogen.sh

GENERIC_CODE = src/lattice.cpp \
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp \
	src/Functional.cpp src/IdealGas.cpp src/ChemicalPotential.cpp \
	src/Minimizer.cpp src/Downhill.cpp src/ExternalPotential.cpp \
	src/QuadraticLineMinimizer.cpp src/SteepestDescent.cpp src/GaussianPolynomial.cpp \
	src/EffectivePotentialToDensity.cpp

bin_PROGRAMS = deft
deft_SOURCES = src/deft.cpp $(GENERIC_CODE)

check_PROGRAMS = \
	tests/simple-liquid.test \
	tests/fftinverse.test \
	tests/ideal-gas.test \
	tests/ideal-gas-at-pressure.test \
	tests/ideal-gas-optimization.test \
	tests/ideal-gas-external-potential.test \
	tests/eps.test

TESTS = tests/null.sh $(check_PROGRAMS)
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

ifeq ($(shell which ikiwiki),)
IKIWIKI=echo "** ikiwiki not found" >&2 ; echo ikiwiki
else
IKIWIKI=ikiwiki
endif

html:
	$(IKIWIKI) `pwd` html -v --wikiname Deft --plugin=goodstuff \
		--exclude=html --exclude=Makefile.am

local-clean:
	rm -rf .ikiwiki html
