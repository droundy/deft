cmake_minimum_required (VERSION 2.6)
project (deft)
include(cmake/flags.cmake)

#add_subdirectory(src/haskell)

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
  src/gSigmaA_by_handFast.cpp src/gSigmaA_automagicFast.cpp
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
include_directories(${CMAKE_CURRENT_BINARY_DIR}/tests)

add_executable(monte-carlo src/Monte-Carlo/monte-carlo.cpp src/utilities.cpp)
add_executable(pair-monte-carlo src/Monte-Carlo/pair-monte-carlo.cpp src/utilities.cpp)


include(cmake/papers.cmake)

# the following are pretty quick programs that we automatically rerun
# on building papers, to ensure our figures are up-to-date.
add_papers_dat_using(deftallcontact
  water-saft/figs/surface-tension
  water-saft/figs/equation-of-state
  water-saft/figs/pressure-with-isotherms

  hughes-saft/figs/surface-tension
  hughes-saft/figs/equation-of-state
  hughes-saft/figs/single-rod-in-water-low-res
  hughes-saft/figs/pressure-with-isotherms

  contact/figs/free-energy
  contact/figs/walls

  pair-correlation/figs/walls
)

# the following are slow programs that we do not automatically rerun
# on building papers
add_papers_mkdat_using(deftallcontact
  water-saft/figs/rods-in-water
  water-saft/figs/single-rod
  water-saft/figs/hughes-single-rod
  water-saft/figs/sphere

  hughes-saft/figs/rods-in-water
  hughes-saft/figs/four-rods-in-water
  hughes-saft/figs/single-rod-in-water-high-res
  hughes-saft/figs/single-rod
  hughes-saft/figs/sphere

  contact/figs/gHS-vs-n
  contact/figs/sphere
  contact/figs/test-particle-wall
  contact/figs/inner-sphere

  pair-correlation/figs/sphere-with-wall

  fuzzy-fmt/figs/walls
)

add_custom_target(
  papers make
  WORKING_DIRECTORY papers
  DEPENDS ${PAPER_DATS}
)

doc_pdf(
  doc/Association.pdf
  doc/WhiteBear.pdf
  doc/TensorWhiteBear.pdf
	doc/WhiteBearMarkII.pdf
  doc/Dispersion.pdf
  doc/SaftFluid.pdf
	doc/SimpDispersion.pdf
  doc/EntropySaftFluid.pdf
	doc/CorrelationS.pdf
  doc/CorrelationA.pdf
  doc/CorrelationAm2.pdf
	doc/GradDispersion.pdf
  doc/JoinedGradDispersion.pdf
  doc/SimpGradDispersion.pdf
)

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
add_simple_tests_for (deftallcontact functional-of-double) # not memory, since that requires some data files
add_simple_tests_for (deftnew new-code)

add_executable(tests/convolve.test tests/convolve.cpp
  src/TensorDensityXXFast.cpp src/VectorDensityXFast.cpp src/n2DensityFast.cpp)
target_link_libraries(tests/convolve.test deftgeneric)
add_test(convolve tests/convolve.test)

# Rules for generated-code.test, which requires some rules to generate
# test code.
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/tests/generated-haskell/nice-sum.h
  COMMAND test -d tests/generated-haskell || mkdir -p tests/generated-haskell
  COMMAND rm -f tests/generated-haskell/*
  COMMAND haskell/test.exe codegen
  DEPENDS haskell/functionals.exe
)

add_executable(tests/generated-code.test tests/generated-code.cpp
  src/TensorWhiteBearFast.cpp tests/generated-haskell/nice-sum.h)
target_link_libraries(tests/generated-code.test defthaskell)
add_test(generated-code tests/generated-code.test)

