EXTRA_DIST = autogen.sh

TESTS = tests/null.sh fftInverseTest epsTest idealGasTest idealGasPressureTest \
	idealGasOptimization

GENERIC_CODE = src/lattice.cpp \
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp \
	src/Functional.cpp src/IdealGas.cpp src/ChemicalPotential.cpp \
	src/Minimizer.cpp src/Downhill.cpp

bin_PROGRAMS = deft
deft_SOURCES = src/deft.cpp $(GENERIC_CODE)

noinst_PROGRAMS = fftInverseTest epsTest idealGasTest idealGasPressureTest \
	idealGasOptimization
fftInverseTest_SOURCES = tests/fftinverse.cpp $(GENERIC_CODE)
epsTest_SOURCES = tests/eps-test.cpp $(GENERIC_CODE)
idealGasTest_SOURCES = tests/ideal-gas.cpp $(GENERIC_CODE)
idealGasPressureTest_SOURCES = tests/ideal-gas-at-pressure.cpp $(GENERIC_CODE)
idealGasOptimization_SOURCES = tests/ideal-gas-optimization.cpp $(GENERIC_CODE)

#noinst_PROGRAMS = optest
#optest_SOURCES = optest.cpp $(GENERIC_CODE)

#BUILT_SOURCES = patchlevel.h

#patchlevel.h: update_patchlevel.sh @DARCS_INVENTORY@ $(DFT_SOURCES) $(CMND_SOURCES)
#	sh update_patchlevel.sh
