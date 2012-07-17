{- This is the white bear version of the fundamental measure theory
functional for the excess free energy of the hard sphere fluid. -}

module WhiteBear
       ( whitebear, correlation_S_WB, correlation_A_WB,
         phi1, phi2, phi3 )
       where

import Expression
import FMT ( n, n0, n1, n2, n2p, n3,
             n1x, n1y, n1z, n2x, n2y, n2z, n2px, n2py, n2pz,
             shell, shell_diam, xshell, yshell, zshell,
             shellPrime, xshellPrime, yshellPrime, zshellPrime,
             rad,
             sqr_n2v, n1v_dot_n2v )

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

phitot :: Expression RealSpace
phitot = var "phitot" "\\Phi" $ phi1 + phi2 + phi3

dwbdn3, dwbdn0, dwbdn2, dwbdn1 :: Expression RealSpace
dwbdn3 =  d (var "phitot" "\\Phi" $ phi1 + phi2 + phi3)-- d phi1 + d phi2 + d phi3
  where d = derive n3 1

dwbdn2 = d (var "phitot" "\\Phi" $ phi1 + phi2 + phi3)-- d phi1 + d phi2 + d phi3
  where d = derive n2 1

dwbdn1 = derive n1 1 (var "phitot" "\\Phi" $ phi1 + phi2 + phi3)

dwbdn0 = derive n0 1 (var "phitot" "\\Phi" $ phi1 + phi2 + phi3)

correlation_S_WB :: Expression RealSpace
correlation_S_WB = var "correlation_S_WB" "g_{\\sigma}^{S}" g_S_WB
    where g_S_WB = dAdR/(kT * n0**2 * 4*pi* (2*rad)**2)
          dAdR = var "dAdR" "\\frac{dA}{dR}" $
                 kT*(dwbdn3*n2
                     - dwbdn2*n2p
                     - dwbdn1*( n2p/(4*pi*rad) + n0)
                     - dwbdn0 * ( n2p/(4*pi*rad**2) + 2*n0/rad )
                     - (derive n2x 1 phitot * n2px + derive n1x 1 phitot * ( n2px + n2x/rad)/(4*pi*rad))
                     - (derive n2y 1 phitot * n2py + derive n1y 1 phitot * ( n2py + n2y/rad)/(4*pi*rad))
                     - (derive n2z 1 phitot * n2pz + derive n1z 1 phitot * ( n2pz + n2z/rad)/(4*pi*rad)))

correlation_A_WB :: Expression RealSpace
correlation_A_WB = var "correlation_A_WB" "g_{\\sigma}^{A}" g_A_WB
    where g_A_WB = dAdR/(kT * n*shell_diam n )
          dAdR = var "dAdR" "\\frac{dA}{dR}" $
                 kT*n*( shell dwbdn3
                        - shellPrime dwbdn2
                        - ( shellPrime dwbdn1 + (shell dwbdn1)/rad ) / (4*pi*rad)
                        - ( shellPrime dwbdn0 + 2*(shell dwbdn0)/rad ) / (4*pi*rad**2)
                        - ( xshellPrime (derive n2x 1 phitot)  +
                            ( xshellPrime (derive n1x 1 phitot) + (xshell (derive n1x 1 phitot))/rad )/(4*pi*rad) )
                        - ( yshellPrime (derive n2y 1 phitot)  +
                            ( yshellPrime (derive n1y 1 phitot) + (yshell (derive n1y 1 phitot))/rad )/(4*pi*rad) )
                        - ( zshellPrime (derive n2z 1 phitot)  +
                            ( zshellPrime (derive n1z 1 phitot) + (zshell (derive n1z 1 phitot))/rad )/(4*pi*rad) ))
