{- This is a simple logarithmic functional for testing purposes. -}

module LogN0 ( log_n0 ) where

import Expression

import FMT ( n0 )

log_n0 :: Expression Scalar
log_n0 = integrate $ log n0
