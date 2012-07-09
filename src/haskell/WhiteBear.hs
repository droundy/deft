{- This is the white bear version of the fundamental measure theory
functional for the excess free energy of the hard sphere fluid. -}

module WhiteBear
       ( whitebear,
         phi1, phi2, phi3 )
       where

import Expression
import FMT ( n0, n1, n2, n3, sqr_n2v, n1v_dot_n2v )

kT :: Type a => Expression a
kT = s_tex "kT" "kT"

vectorThirdTerm :: Expression RealSpace
vectorThirdTerm = n2*(n2**2 - 3*sqr_n2v)

phi1, phi2, phi3 :: Expression RealSpace
phi1 = var "phi1" "\\Phi_1" $ -n0*log(1-n3)
phi2 = var "phi2" "\\Phi_2" $ (n2*n1 - n1v_dot_n2v)/(1-n3)
phi3 = var "phi3" "\\Phi_3" $ (n3 + (1-n3)**2*log(1-n3))/(36*pi * n3**2 * (1-n3)**2)*vectorThirdTerm


whitebear :: Expression Scalar
whitebear = var "whitebear" "F_{HS}" $ integrate (kT*(phi1+phi2+phi3))
