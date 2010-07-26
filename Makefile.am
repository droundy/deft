EXTRA_DIST = autogen.sh

TESTS = tests/null.sh fftInverseTest

GENERIC_CODE = src/lattice.cpp \
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp

bin_PROGRAMS = deft
deft_SOURCES = src/deft.cpp $(GENERIC_CODE)

noinst_PROGRAMS = fftInverseTest
fftInverseTest_SOURCES = tests/fftinverse.cpp $(GENERIC_CODE)

#noinst_PROGRAMS = optest
#optest_SOURCES = optest.cpp $(GENERIC_CODE)

#BUILT_SOURCES = patchlevel.h

#patchlevel.h: update_patchlevel.sh @DARCS_INVENTORY@ $(DFT_SOURCES) $(CMND_SOURCES)
#	sh update_patchlevel.sh
