{- This is a simple quadratic functional for testing purposes. -}

module Quadratic ( quadratic ) where

import Expression

quadratic :: Expression Scalar
quadratic = integrate $ 0.5*konst*(x-meanval)**2
  where x = r_var "x"
        meanval = s_var "meanval"
        konst = s_var "konst"
