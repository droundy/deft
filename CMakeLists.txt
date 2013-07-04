cmake_minimum_required (VERSION 2.6)
project (deft)
include(cmake/flags.cmake)

include(cmake/haskellcode.cmake)
haskell_functionals(
  src/HardSpheresNoTensor2Fast.cpp
	src/TensorWhiteBearFast.cpp
  src/WhiteBearMarkIIFast.cpp
	src/YuWuCorrelationFast.cpp
  src/SaftFluid2Fast.cpp
	src/EntropySaftFluid2Fast.cpp
  src/CorrelationGrossCorrectFast.cpp
	src/gSigmaSm2Fast.cpp
  src/gSigmaAm2Fast.cpp
	src/gSigmaSFast.cpp
  src/gSigmaAFast.cpp
	src/SoftFluidFast.cpp
  src/HardFluidFast.cpp
	src/WaterXFast.cpp
  src/HughesXFast.cpp
  # the following are only needed for tests
  src/TensorDensityXXFast.cpp src/VectorDensityXFast.cpp src/n2DensityFast.cpp
)

new_haskell_functionals(
  src/new/WhiteBearFast.cpp
)
add_library(deftnew STATIC
  src/new/Minimizer.cpp
)
target_link_libraries(deftnew deftnewhaskell)

add_library(deftgeneric STATIC
  src/lattice.cpp src/utilities.cpp
	src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp
	src/IdealGas.cpp src/ChemicalPotential.cpp
	src/HardSpheres.cpp src/ExternalPotential.cpp
	src/Functional.cpp
	src/Gaussian.cpp src/Pow.cpp
  src/EffectivePotentialToDensity.cpp
	src/equation-of-state.cpp src/water-constants.cpp
	src/compute-surface-tension.cpp
	src/Minimizer.cpp src/Downhill.cpp
	src/Precision.cpp src/ConjugateGradient.cpp
	src/WaterSaftFast.cpp
	src/QuadraticLineMinimizer.cpp src/SteepestDescent.cpp)
target_link_libraries(deftgeneric fftw3) # need ffw3!

add_library(deftcontact STATIC src/ContactDensity.cpp)
target_link_libraries(deftcontact deftgeneric)

add_library(deftallcontact STATIC src/ContactDensity.cpp)
target_link_libraries(deftallcontact defthaskell)

include_directories(include)
include_directories(src)

add_executable(monte-carlo src/Monte-Carlo/monte-carlo.cpp src/utilities.cpp)
add_executable(pair-monte-carlo src/Monte-Carlo/pair-monte-carlo.cpp src/utilities.cpp)

include(cmake/tests.cmake)

add_test(run-pair-monte-carlo pair-monte-carlo 10 1000000
  0.01 /tmp/foo /tmp/dafoo periodxy 20 wallz 20 flatdiv)
add_test(run-monte-carlo monte-carlo 10 100000 0.01 /tmp/test.out)

add_simple_tests_for (deftgeneric surface-tension functional-arithmetic)
add_simple_tests_for (defthaskell
  saft eos eps fftinverse ideal-gas precision
  print-iter convolve-finite-difference
)
add_simple_tests_for (deftnewhaskell new-hard-spheres)
add_simple_tests_for (deftallcontact memory functional-of-double)
add_simple_tests_for (deftnew new-code)

add_executable(tests/convolve.test tests/convolve.cpp
  src/TensorDensityXXFast.cpp src/VectorDensityXFast.cpp src/n2DensityFast.cpp)
target_link_libraries(tests/convolve.test deftgeneric)
add_test(convolve tests/convolve.test)

# Rules for generated-code.test, which requires some rules to generate
# test code.
add_custom_command(
  OUTPUT tests/generated-haskell/nice-sum.h
  COMMAND test -d tests/generated-haskell || mkdir tests/generated-haskell
  COMMAND rm -f tests/generated-haskell/*
  COMMAND src/haskell/test.exe codegen
  DEPENDS src/haskell/functionals.exe
)
add_executable(tests/generated-code.test tests/generated-code.cpp
  src/TensorWhiteBearFast.cpp tests/generated-haskell/nice-sum.h)
target_link_libraries(tests/generated-code.test deftgeneric)
add_test(generated-code tests/generated-code.test)

