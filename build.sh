#!/bin/sh

set -ev

(cd src/haskell && python2 create_generators.py)

(cd src/haskell && ghc -O2 -c LatexDouble.hs)

(cd src/haskell && ghc -O2 -c Expression.hs)

(cd src/haskell && ghc -O2 -c C.hs)

(cd src/haskell && ghc -O2 -c Statement.hs)

(cd src/haskell && ghc -O2 -c Optimize.hs)

(cd src/haskell && ghc -O2 -c CodeGen.hs)

(cd src/haskell && ghc -O2 -c IdealGas.hs)

(cd src/haskell && ghc -O2 -c FMT.hs)

(cd src/haskell && ghc -O2 -c WhiteBear.hs)

(cd src/haskell && ghc -O2 -c WaterSaft.hs)

(cd src/haskell && ghc -O2 -c HughesSaft.hs)

(cd src/haskell && ghc -O2 -c ExternalLennardJones.hs)

(cd src/haskell && ghc -O2 -c ExternalPotentialTest.hs)

(cd src/haskell && ghc -O2 -c SFMT.hs)

(cd src/haskell && ghc -O2 -c Rosenfeld.hs)

(cd src/haskell && ghc -O2 -c functionals.hs)

(cd src/haskell && ghc -O2 -c Latex.hs)

(cd src/haskell && ghc -O2 -c latex-functionals.hs)

(cd src/haskell && ghc -O2 -c LogN0.hs)

(cd src/haskell && ghc -O2 -c NewCode.hs)

(cd src/haskell && ghc -O2 -c QuadraticGaussian.hs)

(cd src/haskell && ghc -O2 -c Quadratic.hs)

(cd src/haskell && ghc -O2 -c QuadraticN0.hs)

(cd src/haskell && ghc -O2 -c test.hs)

(cd src/haskell && ghc -O2 --make -o generate_ExternalPotentialTest.exe generate_ExternalPotentialTest.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_ExternalPotentialTest.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/ExternalPotentialTestFast.o src/new/ExternalPotentialTestFast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_HomogeneousSFMTFluid.exe generate_HomogeneousSFMTFluid.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_HomogeneousSFMTFluid.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/HomogeneousSFMTFluidFast.o src/new/HomogeneousSFMTFluidFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/ExternalPotential.o src/ExternalPotential.cpp)

(cd src/haskell && ghc -O2 --make -o generate_HomogeneousWaterSaftByHand.exe generate_HomogeneousWaterSaftByHand.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_HomogeneousWaterSaftByHand.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/HomogeneousWaterSaftByHandFast.o src/new/HomogeneousWaterSaftByHandFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o papers/fuzzy-fmt/figs/radial-wca.o papers/fuzzy-fmt/figs/radial-wca.cpp)

(cd src/haskell && ghc -O2 --make -o generate_HomogeneousWaterSaft.exe generate_HomogeneousWaterSaft.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_HomogeneousWaterSaft.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/HomogeneousWaterSaftFast.o src/new/HomogeneousWaterSaftFast.cpp)

(cd src/haskell && ghc -O2 -package containers -package filepath -package directory -o functionals.exe SFMT.o WhiteBear.o Rosenfeld.o WaterSaft.o HughesSaft.o LatexDouble.o FMT.o functionals.o IdealGas.o CodeGen.o Statement.o Expression.o Optimize.o)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/gSigmaA_automagicFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/gSigmaA_automagicFast.o src/gSigmaA_automagicFast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_HomogeneousWhiteBear.exe generate_HomogeneousWhiteBear.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_HomogeneousWhiteBear.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/HomogeneousWhiteBearFast.o src/new/HomogeneousWhiteBearFast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_QuadraticGaussian.exe generate_QuadraticGaussian.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_QuadraticGaussian.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/QuadraticGaussianFast.o src/new/QuadraticGaussianFast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_LogN0.exe generate_LogN0.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_LogN0.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/LogN0Fast.o src/new/LogN0Fast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_SFMTFluid.exe generate_SFMTFluid.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_SFMTFluid.exe)

(cd src/haskell && ghc -O2 --make -o generate_SFMTFluidVeff.exe generate_SFMTFluidVeff.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_SFMTFluidVeff.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o papers/fuzzy-fmt/figs/new-melting.o papers/fuzzy-fmt/figs/new-melting.cpp)

(cd src/haskell && ghc -O2 --make -o generate_Phi1.exe generate_Phi1.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_Phi1.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/Phi1Fast.o src/new/Phi1Fast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_SPhi1.exe generate_SPhi1.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_SPhi1.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/SPhi1Fast.o src/new/SPhi1Fast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/TensorWhiteBearFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/TensorWhiteBearFast.o src/TensorWhiteBearFast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_Phi2.exe generate_Phi2.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_Phi2.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/Phi2Fast.o src/new/Phi2Fast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_Phi3.exe generate_Phi3.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_Phi3.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/Phi3Fast.o src/new/Phi3Fast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/radial-distribution-monte-carlo.o src/Monte-Carlo/radial-distribution-monte-carlo.cpp)

(cd src/haskell && ghc -O2 --make -o generate_Quadratic.exe generate_Quadratic.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_Quadratic.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/QuadraticFast.o src/new/QuadraticFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Functional.o src/Functional.cpp)

(cd src/haskell && ghc -O2 --make -o generate_QuadraticN0.exe generate_QuadraticN0.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_QuadraticN0.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/QuadraticN0Fast.o src/new/QuadraticN0Fast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/SFMTFluidFast.o src/new/SFMTFluidFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/ChemicalPotential.o src/ChemicalPotential.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/SFMTFluidVeffFast.o src/new/SFMTFluidVeffFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/gSigmaA_by_handFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/gSigmaA_by_handFast.o src/gSigmaA_by_handFast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_SPhi2.exe generate_SPhi2.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_SPhi2.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/SPhi2Fast.o src/new/SPhi2Fast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/VectorDensityXFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/VectorDensityXFast.o src/VectorDensityXFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Faddeeva.o src/Faddeeva.cpp)

(cd src/haskell && ghc -O2 --make -o generate_SPhi3.exe generate_SPhi3.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_SPhi3.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/SPhi3Fast.o src/new/SPhi3Fast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_WaterSaftByHand.exe generate_WaterSaftByHand.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_WaterSaftByHand.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/WaterSaftByHandFast.o src/new/WaterSaftByHandFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/free-energy-monte-carlo.o src/Monte-Carlo/free-energy-monte-carlo.cpp)

(cd src/haskell && ghc -O2 --make -o generate_WaterSaft.exe generate_WaterSaft.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_WaterSaft.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/WaterSaftFast.o src/new/WaterSaftFast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_WhiteBearFluid.exe generate_WhiteBearFluid.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_WhiteBearFluid.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/WhiteBearFluidFast.o src/new/WhiteBearFluidFast.cpp)

(cd src/haskell && ghc -O2 --make -o generate_WhiteBear.exe generate_WhiteBear.hs)

(cd src/haskell && cd ../.. && src/haskell/generate_WhiteBear.exe)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/WhiteBearFast.o src/new/WhiteBearFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/HughesXFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/HughesXFast.o src/HughesXFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/polyhedra.o src/Monte-Carlo/polyhedra.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/HardSpheresNoTensor2Fast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/HardSpheresNoTensor2Fast.o src/HardSpheresNoTensor2Fast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/WhiteBearMarkIIFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/WhiteBearMarkIIFast.o src/WhiteBearMarkIIFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/TensorDensityXXFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/TensorDensityXXFast.o src/TensorDensityXXFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/SaftFluid2Fast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/SaftFluid2Fast.o src/SaftFluid2Fast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/WaterSaftFast.o src/WaterSaftFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/n2DensityFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/n2DensityFast.o src/n2DensityFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/YuWuCorrelationFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/YuWuCorrelationFast.o src/YuWuCorrelationFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/gSigmaAFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/gSigmaAFast.o src/gSigmaAFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Minimizer.o src/Minimizer.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/EntropySaftFluid2Fast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/EntropySaftFluid2Fast.o src/EntropySaftFluid2Fast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/CorrelationGrossCorrectFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/CorrelationGrossCorrectFast.o src/CorrelationGrossCorrectFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/gSigmaSm2Fast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/gSigmaSm2Fast.o src/gSigmaSm2Fast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/gSigmaAm2Fast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/gSigmaAm2Fast.o src/gSigmaAm2Fast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/gSigmaSFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/gSigmaSFast.o src/gSigmaSFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/HughesHBFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/HughesHBFast.o src/HughesHBFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/SoftFluidFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/SoftFluidFast.o src/SoftFluidFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/square-well.o src/Monte-Carlo/square-well.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/HardFluidFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/HardFluidFast.o src/HardFluidFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/HardRosenfeldFluidFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/HardRosenfeldFluidFast.o src/HardRosenfeldFluidFast.cpp)

(cd src/haskell && cd ../.. && src/haskell/functionals.exe src/WaterXFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/WaterXFast.o src/WaterXFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/ReciprocalGrid.o src/ReciprocalGrid.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/HardSpheres.o src/HardSpheres.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/lattice.o src/lattice.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/utilities.o src/utilities.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/WaterSaft_by_handFast.o src/WaterSaft_by_handFast.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/IdealGas.o src/IdealGas.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/ContactDensity.o src/ContactDensity.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/NewFunctional.o src/new/NewFunctional.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/ConjugateGradient.o src/ConjugateGradient.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Grid.o src/Grid.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Gaussian.o src/Gaussian.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/SteepestDescent.o src/SteepestDescent.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/GridDescription.o src/GridDescription.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/vector3d.o src/vector3d.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Precision.o src/Precision.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Pow.o src/Pow.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/equation-of-state.o src/equation-of-state.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/new/Minimize.o src/new/Minimize.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o papers/fuzzy-fmt/figs/homogeneous.o papers/fuzzy-fmt/figs/homogeneous.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/compute-surface-tension.o src/compute-surface-tension.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/EffectivePotentialToDensity.o src/EffectivePotentialToDensity.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/QuadraticLineMinimizer.o src/QuadraticLineMinimizer.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Downhill.o src/Downhill.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/water-constants.o src/water-constants.cpp)

(g++ -lpopt -lprofiler -lfftw3 -o papers/fuzzy-fmt/figs/homogeneous.mkdat src/ReciprocalGrid.o src/HardSpheres.o src/lattice.o src/ExternalPotential.o src/utilities.o src/WaterSaft_by_handFast.o src/IdealGas.o src/ContactDensity.o src/new/NewFunctional.o src/ConjugateGradient.o src/Grid.o src/Gaussian.o src/WaterSaftFast.o src/SteepestDescent.o src/GridDescription.o src/vector3d.o src/Precision.o src/Pow.o src/Functional.o src/equation-of-state.o src/HardFluidFast.o src/new/Minimize.o papers/fuzzy-fmt/figs/homogeneous.o src/compute-surface-tension.o src/Faddeeva.o src/ChemicalPotential.o src/SoftFluidFast.o src/EffectivePotentialToDensity.o src/QuadraticLineMinimizer.o src/Downhill.o src/water-constants.o src/Minimizer.o)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o papers/fuzzy-fmt/figs/new-walls.o papers/fuzzy-fmt/figs/new-walls.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/monte-carlo.o src/Monte-Carlo/monte-carlo.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/soft-monte-carlo.o src/Monte-Carlo/soft-monte-carlo.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/pair-monte-carlo.o src/Monte-Carlo/pair-monte-carlo.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/triplet-monte-carlo.o src/Monte-Carlo/triplet-monte-carlo.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/polyhedra-monte-carlo.o src/Monte-Carlo/polyhedra-monte-carlo.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/polyhedra-talk.o src/Monte-Carlo/polyhedra-talk.cpp)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o src/Monte-Carlo/square-well-monte-carlo.o src/Monte-Carlo/square-well-monte-carlo.cpp)

(g++ -lpopt -lprofiler -lfftw3 -o monte-carlo src/utilities.o src/Monte-Carlo/monte-carlo.o)

(g++ -lpopt -lprofiler -lfftw3 -o soft-monte-carlo src/utilities.o src/Monte-Carlo/soft-monte-carlo.o)

(g++ -lpopt -lprofiler -lfftw3 -o pair-monte-carlo src/Monte-Carlo/pair-monte-carlo.o src/utilities.o)

(g++ -lpopt -lprofiler -lfftw3 -o triplet-monte-carlo src/utilities.o src/Monte-Carlo/triplet-monte-carlo.o)

(g++ -lpopt -lprofiler -lfftw3 -o polyhedra-monte-carlo src/vector3d.o src/Monte-Carlo/polyhedra-monte-carlo.o src/Monte-Carlo/polyhedra.o)

(g++ -lpopt -lprofiler -lfftw3 -o polyhedra-talk src/vector3d.o src/Monte-Carlo/polyhedra-talk.o src/Monte-Carlo/polyhedra.o)

(g++ -lpopt -lprofiler -lfftw3 -o square-well-monte-carlo src/Monte-Carlo/square-well.o src/vector3d.o src/Monte-Carlo/square-well-monte-carlo.o)

(g++ -lpopt -lprofiler -lfftw3 -o radial-distribution-monte-carlo src/Monte-Carlo/radial-distribution-monte-carlo.o src/utilities.o)

(g++ -lpopt -lprofiler -lfftw3 -o free-energy-monte-carlo src/Monte-Carlo/square-well.o src/Monte-Carlo/free-energy-monte-carlo.o src/vector3d.o)

(g++ -lpopt -lprofiler -lfftw3 -o papers/fuzzy-fmt/figs/radial-wca.mkdat src/ReciprocalGrid.o src/HardSpheres.o src/lattice.o src/ExternalPotential.o src/utilities.o src/WaterSaft_by_handFast.o src/IdealGas.o papers/fuzzy-fmt/figs/radial-wca.o src/ContactDensity.o src/new/NewFunctional.o src/ConjugateGradient.o src/Grid.o src/Gaussian.o src/WaterSaftFast.o src/SteepestDescent.o src/GridDescription.o src/vector3d.o src/Precision.o src/Pow.o src/Functional.o src/equation-of-state.o src/new/Minimize.o src/Faddeeva.o src/compute-surface-tension.o src/ChemicalPotential.o src/SoftFluidFast.o src/EffectivePotentialToDensity.o src/QuadraticLineMinimizer.o src/Downhill.o src/water-constants.o src/Minimizer.o)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o papers/fuzzy-fmt/figs/soft-wall.o papers/fuzzy-fmt/figs/soft-wall.cpp)

(g++ -lpopt -lprofiler -lfftw3 -o papers/fuzzy-fmt/figs/soft-wall.mkdat papers/fuzzy-fmt/figs/soft-wall.o src/ReciprocalGrid.o src/HardSpheres.o src/lattice.o src/ExternalPotential.o src/utilities.o src/WaterSaft_by_handFast.o src/IdealGas.o src/ContactDensity.o src/new/NewFunctional.o src/ConjugateGradient.o src/Grid.o src/Gaussian.o src/WaterSaftFast.o src/SteepestDescent.o src/GridDescription.o src/vector3d.o src/Precision.o src/Pow.o src/Functional.o src/equation-of-state.o src/new/Minimize.o src/Faddeeva.o src/compute-surface-tension.o src/ChemicalPotential.o src/SoftFluidFast.o src/EffectivePotentialToDensity.o src/QuadraticLineMinimizer.o src/Downhill.o src/water-constants.o src/Minimizer.o)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o papers/fuzzy-fmt/figs/new-soft-wall.o papers/fuzzy-fmt/figs/new-soft-wall.cpp)

(g++ -lpopt -lprofiler -lfftw3 -o papers/fuzzy-fmt/figs/new-soft-wall.mkdat src/ReciprocalGrid.o src/HardSpheres.o src/lattice.o src/ExternalPotential.o src/new/SFMTFluidVeffFast.o src/WaterSaft_by_handFast.o src/IdealGas.o src/equation-of-state.o src/ContactDensity.o src/new/NewFunctional.o papers/fuzzy-fmt/figs/new-soft-wall.o src/ConjugateGradient.o src/Grid.o src/Gaussian.o src/WaterSaftFast.o src/SteepestDescent.o src/GridDescription.o src/vector3d.o src/Precision.o src/Pow.o src/new/SFMTFluidFast.o src/Functional.o src/utilities.o src/new/Minimize.o src/Faddeeva.o src/compute-surface-tension.o src/ChemicalPotential.o src/EffectivePotentialToDensity.o src/QuadraticLineMinimizer.o src/Downhill.o src/new/HomogeneousSFMTFluidFast.o src/water-constants.o src/Minimizer.o)

(g++ -lpopt -lprofiler -lfftw3 -o papers/fuzzy-fmt/figs/new-melting.mkdat src/ReciprocalGrid.o src/HardSpheres.o src/lattice.o src/ExternalPotential.o src/new/SFMTFluidVeffFast.o src/WaterSaft_by_handFast.o src/IdealGas.o src/equation-of-state.o src/ContactDensity.o src/new/NewFunctional.o src/ConjugateGradient.o src/Grid.o src/Gaussian.o src/WaterSaftFast.o src/SteepestDescent.o src/GridDescription.o src/vector3d.o src/Precision.o src/Pow.o src/new/SFMTFluidFast.o papers/fuzzy-fmt/figs/new-melting.o src/Functional.o src/utilities.o src/new/Minimize.o src/Faddeeva.o src/compute-surface-tension.o src/ChemicalPotential.o src/EffectivePotentialToDensity.o src/QuadraticLineMinimizer.o src/Downhill.o src/new/HomogeneousSFMTFluidFast.o src/water-constants.o src/Minimizer.o)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o papers/fuzzy-fmt/figs/walls.o papers/fuzzy-fmt/figs/walls.cpp)

(g++ -lpopt -lprofiler -lfftw3 -o papers/fuzzy-fmt/figs/walls.mkdat src/ReciprocalGrid.o src/HardSpheres.o src/lattice.o papers/fuzzy-fmt/figs/walls.o src/ExternalPotential.o src/utilities.o src/WaterSaft_by_handFast.o src/IdealGas.o src/ContactDensity.o src/new/NewFunctional.o src/ConjugateGradient.o src/Grid.o src/Gaussian.o src/WaterSaftFast.o src/SteepestDescent.o src/GridDescription.o src/vector3d.o src/HardRosenfeldFluidFast.o src/Precision.o src/Pow.o src/Functional.o src/equation-of-state.o src/new/Minimize.o src/Faddeeva.o src/compute-surface-tension.o src/ChemicalPotential.o src/SoftFluidFast.o src/EffectivePotentialToDensity.o src/QuadraticLineMinimizer.o src/Downhill.o src/water-constants.o src/Minimizer.o)

(g++ -lpopt -lprofiler -lfftw3 -o papers/fuzzy-fmt/figs/new-walls.mkdat src/ReciprocalGrid.o src/HardSpheres.o src/lattice.o src/ExternalPotential.o src/new/SFMTFluidVeffFast.o src/WaterSaft_by_handFast.o src/IdealGas.o src/equation-of-state.o src/ContactDensity.o src/new/NewFunctional.o src/ConjugateGradient.o src/Grid.o src/Gaussian.o src/WaterSaftFast.o src/SteepestDescent.o src/GridDescription.o src/vector3d.o src/Precision.o src/Pow.o src/new/SFMTFluidFast.o src/Functional.o src/utilities.o src/new/Minimize.o papers/fuzzy-fmt/figs/new-walls.o src/Faddeeva.o src/compute-surface-tension.o src/ChemicalPotential.o src/EffectivePotentialToDensity.o src/QuadraticLineMinimizer.o src/Downhill.o src/new/HomogeneousSFMTFluidFast.o src/water-constants.o src/Minimizer.o)

(g++ -Wall -Werror -O3 -std=c++11 -Isrc -Iinclude -g -c -o papers/fuzzy-fmt/figs/soft-sphere.o papers/fuzzy-fmt/figs/soft-sphere.cpp)

(g++ -lpopt -lprofiler -lfftw3 -o papers/fuzzy-fmt/figs/soft-sphere.mkdat src/ReciprocalGrid.o src/HardSpheres.o src/lattice.o src/ExternalPotential.o src/utilities.o src/WaterSaft_by_handFast.o src/IdealGas.o src/ContactDensity.o src/new/NewFunctional.o src/ConjugateGradient.o src/Grid.o src/Gaussian.o src/WaterSaftFast.o src/SteepestDescent.o src/GridDescription.o src/vector3d.o src/Precision.o src/Pow.o src/Functional.o src/equation-of-state.o src/new/Minimize.o src/Faddeeva.o src/compute-surface-tension.o src/ChemicalPotential.o src/SoftFluidFast.o src/EffectivePotentialToDensity.o src/QuadraticLineMinimizer.o src/Downhill.o papers/fuzzy-fmt/figs/soft-sphere.o src/water-constants.o src/Minimizer.o)

