{- This is a simple quadratic functional for testing purposes. -}

module QuadraticGaussian ( quadratic_gaussian ) where

import Expression

quadratic_gaussian :: Expression Scalar
quadratic_gaussian = integrate $ 0.5*konst*(gaussian-meanval)**2
  where meanval = s_var "meanval"
        konst = s_var "konst"

gaussian :: Expression RealSpace
gaussian = ifft ( gaussiank * fft x)
  where gaussiank = exp(-(k*sigma)**2/2)
        sigma = s_var "sigma"
        x = r_var "x"
