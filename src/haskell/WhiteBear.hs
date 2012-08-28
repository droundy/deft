{- This is the white bear version of the fundamental measure theory
functional for the excess free energy of the hard sphere fluid. -}

module WhiteBear
       ( whitebear, correlation_S_WB, correlation_A_WB,
         whitebear_m2, correlation_S_WB_m2, correlation_A_WB_m2,
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

correlation_S_WB :: Expression RealSpace
correlation_S_WB = var "correlation_S_WB" "g_{\\sigma}^{S}" $
                   correlation_S_helper phitot

correlation_S_WB_m2 :: Expression RealSpace
correlation_S_WB_m2 = var "correlation_S_WB" "g_{\\sigma,m2}^{S}" $
                      correlation_S_helper phitot_m2

correlation_S_helper :: Expression RealSpace -> Expression RealSpace
correlation_S_helper phit = dAdR/(kT * n0**2 * 4*pi* (2*rad)**2)
    where dAdR = var "dAdR" "\\frac{dA}{dR}" $
                 kT*((d n3)*n2
                     - (d n2)*n2p
                     - (d n1)*( n2p/(4*pi*rad) + n0)
                     - (d n0) * ( n2p/(4*pi*rad**2) + 2*n0/rad )
                     - (d n2x * n2px + d n1x * ( n2px - n2x/rad)/(4*pi*rad))
                     - (d n2y * n2py + d n1y * ( n2py - n2y/rad)/(4*pi*rad))
                     - (d n2z * n2pz + d n1z * ( n2pz - n2z/rad)/(4*pi*rad)))
          d a = derive a 1 phit

correlation_A_WB :: Expression RealSpace
correlation_A_WB = var "correlation_A_WB" "g_{\\sigma}^{A}" $
                   correlation_A_helper phitot

correlation_A_WB_m2 :: Expression RealSpace
correlation_A_WB_m2 = var "correlation_A_WB" "g_{\\sigma,m2}^{A}" $
                      correlation_A_helper phitot_m2

correlation_A_helper :: Expression RealSpace -> Expression RealSpace
correlation_A_helper phit = dAdR/(kT * n*shell_diam n )
    where dAdR = var "dAdR" "\\frac{dA}{dR}" $
                 kT*n*( shell (d n3)
                        - shellPrime (d n2)
                        - ( shellPrime (d n1) + (shell (d n1))/rad ) / (4*pi*rad)
                        - ( shellPrime (d n0) + 2*(shell (d n0))/rad ) / (4*pi*rad**2)
                        + ( xshellPrime (d n2x)  +
                            ( xshellPrime (d n1x) + (xshell (d n1x))/rad )/(4*pi*rad) )
                        + ( yshellPrime (d n2y)  +
                            ( yshellPrime (d n1y) + (yshell (d n1y))/rad )/(4*pi*rad) )
                        + ( zshellPrime (d n2z)  +
                            ( zshellPrime (d n1z) + (zshell (d n1z))/rad )/(4*pi*rad) ))
          d a = derive a 1 phit


whitebear_m2 :: Expression Scalar
whitebear_m2 = var "whitebear_m2" "F_{HS}^{(m2)}" $ integrate (kT*phitot_m2)

phi2_m2, phi3_m2, phitot_m2 :: Expression RealSpace
phitot_m2 = var "phitot_m2" "\\Phi^{(m2)}" $
            phi1 + phi2_m2 + phi3_m2

phi2_m2 = var "Phi2ii" "\\Phi_2^{(m2)}" $
          (1+(1/9)*(n3**2)*phi2II)*(n2*n1 - n1v_dot_n2v)/(1-n3)
            where -- above is from equation 16 of Goos 2006
                  -- below is from equation 17 of Goos 2006
                  phi2II = --var "phi2ii" "\\phi_2^{(m2)}" $
                           (6*n3 - 3*n3**2 + 6*(1-n3)*log(1-n3))/(n3**3)
                    --(1 + n3/2)
phi3_m2 = var "Phi3ii" "\\Phi_3^{(m2)}" $
          (1-(4/9)*n3*phi3II)*(n2**3 - 3*n2*sqr_n2v)/(24*pi*(1-n3)**2)
            where -- above is from equation 16 of Goos 2006
                  -- below is from equation 17 of Goos 2006
                  phi3II = var "phi3ii" "\\phi_3^{(m2)}" $
                           (6*n3 - 9*n3**2 + 6*n3**3 + 6*(1-n3)**2*log(1-n3))/(4*n3**3)
