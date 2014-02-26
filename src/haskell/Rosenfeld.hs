
module Rosenfeld
       ( fmt )
       where

import Expression
import WhiteBear ( kT )
import FMT ( n0, n1, n2, n3, sqr_n2v, n1v_dot_n2v )

phi1, phi2, phi3 :: Expression Scalar
phi1 = var "phi1" "\\Phi_1" $ integrate $ -kT*n0*log(1-n3)
phi2 = var "phi2" "\\Phi_2" $ integrate $ kT*(n2*n1 - n1v_dot_n2v)/(1-n3)
phi3 = var "phi3" "\\Phi_3" $ integrate $ kT*(n2**3/3 - sqr_n2v*n2)/(8*pi*(1-n3)**2)

fmt :: Expression Scalar
fmt = var "fmt" "F_{\\text{hard}}" $ (phi1 + phi2 + phi3)
