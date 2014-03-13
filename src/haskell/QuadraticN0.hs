{- This is a simple quadratic functional for testing purposes. -}

module QuadraticN0 ( quadratic_n0 ) where

import Expression

import FMT ( n0 )

quadratic_n0 :: Expression Scalar
quadratic_n0 = integrate $ 0.5*konst*(n0-meanval)**2
  where meanval = s_var "meanval"
        konst = s_var "konst"
