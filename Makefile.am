EXTRA_DIST = autogen.sh

GENERIC_CODE = src/lattice.cpp \
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp \
	src/Functional.cpp src/IdealGas.cpp src/ChemicalPotential.cpp \
	src/Minimizer.cpp src/Downhill.cpp

bin_PROGRAMS = deft
deft_SOURCES = src/deft.cpp $(GENERIC_CODE)

check_PROGRAMS = \
	tests/fftinverse.test \
	tests/eps.test \
	tests/ideal-gas.test \
	tests/ideal-gas-at-pressure.test \
	tests/ideal-gas-optimization.test
TESTS = tests/null.sh $(check_PROGRAMS)
tests_fftinverse_test_SOURCES = tests/fftinverse.cpp $(GENERIC_CODE)
tests_eps_test_SOURCES = tests/eps.cpp $(GENERIC_CODE)
tests_ideal_gas_test_SOURCES = tests/ideal-gas.cpp $(GENERIC_CODE)
tests_ideal_gas_at_pressure_test_SOURCES = tests/ideal-gas-at-pressure.cpp $(GENERIC_CODE)
tests_ideal_gas_optimization_test_SOURCES = tests/ideal-gas-optimization.cpp $(GENERIC_CODE)

#noinst_PROGRAMS = optest
#optest_SOURCES = optest.cpp $(GENERIC_CODE)

#BUILT_SOURCES = patchlevel.h

#patchlevel.h: update_patchlevel.sh @DARCS_INVENTORY@ $(DFT_SOURCES) $(CMND_SOURCES)
#	sh update_patchlevel.sh
