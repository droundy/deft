#!/usr/bin/python2

for x in ['HardSpheresNoTensor2Fast',
          'TensorWhiteBearFast',
          'WhiteBearMarkIIFast',
          'gSigmaA_by_handFast',
          'gSigmaA_automagicFast',
          'TensorDensityXXFast',
          'n2DensityFast',
          'VectorDensityXFast',
          'YuWuCorrelationFast',
          'SaftFluid2Fast',
          'EntropySaftFluid2Fast',
          'CorrelationGrossCorrectFast',
          'gSigmaSm2Fast',
          'gSigmaAm2Fast',
          'gSigmaSFast',
          'gSigmaAFast',
          'HughesHBFast',
          'SoftFluidFast',
          'HardFluidFast',
          'HardRosenfeldFluidFast',
          'WaterXFast',
          'HughesXFast']:
    print ': haskell/functionals.exe |> cd .. && src/%%f src/%%o |> %s' % (x + '.cpp')
    #print ': %s.cpp |> !cxx |>' % x
