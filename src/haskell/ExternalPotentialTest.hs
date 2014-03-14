{- This is a simple functional for testing purposes. -}

module ExternalPotentialTest ( external_potential ) where

import Expression

external_potential :: Expression Scalar
external_potential = integrate $ v*n
  where v = r_var "V"
        n = r_var "n"
