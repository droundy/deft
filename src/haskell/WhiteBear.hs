{- This is the white bear version of the fundamental measure theory
functional for the excess free energy of the hard sphere fluid. -}

module WhiteBear
       ( whitebear_n, homogeneous_whitebear,
         whitebear_fluid_n, whitebear_fluid_Veff, homogeneous_whitebear_fluid,
         kT, whitebear, gSigmaS, gSigmaA, gSigmaA_automagic, gSigmaA_by_hand,
         tensorThirdTerm, phi3t,
         tensorwhitebear,
         whitebear_m2, gSigmaS_m2, gSigmaA_m2,
         kTphi1, kTphi2, kTphi3,
         correlation_gross, phi1, phi2, phi3, nA )
       where

import Expression
import FMT ( n, n0, n1, n2, n3, n1v, n2v, n2m,
             shell, shell_diam, step_diam, vshelldot,
             shellPrime, vshellPrimedot,
             rad, smear, kR,
             sqr_n2v, n1v_dot_n2v )
import IdealGas ( kT, idealgas )

whitebear_n, homogeneous_whitebear, whitebear_fluid_n, whitebear_fluid_Veff :: Expression Scalar
whitebear_n = substitute n (r_var "n") whitebear
homogeneous_whitebear = makeHomogeneous whitebear_n
whitebear_fluid_n = substitute n (r_var "n") $
                    whitebear + idealgas + integrate (n * (r_var "Vext" - s_var "mu"))
whitebear_fluid_Veff = substitute n ("n" === exp(- r_var "Veff"/kT)) $
                       whitebear + idealgas + integrate (n * (r_var "Vext" - s_var "mu"))

homogeneous_whitebear_fluid :: Expression Scalar
homogeneous_whitebear_fluid = makeHomogeneous $ substitute n (r_var "n") $
                              whitebear + idealgas - integrate (n * s_var "mu")


vectorThirdTerm :: Expression RealSpace
vectorThirdTerm = n2*(n2**2 - 3*sqr_n2v)

tensorThirdTerm :: Expression RealSpace
tensorThirdTerm = n2**3 - 3*n2*sqr_n2v + 9*((n2v .*. n2m `dot` n2v) - tracetensorcube n2m/2)

phi1, phi2, phi3, phi3t :: Expression RealSpace
phi1 = var "phi1" "\\Phi_1" $ -n0*log(1-n3)
phi2 = var "phi2" "\\Phi_2" $ (n2*n1 - n1v_dot_n2v)/(1-n3)
phi3 = var "phi3" "\\Phi_3" $ (n3 + (1-n3)**2*log(1-n3))/(36*pi * n3**2 * (1-n3)**2)*vectorThirdTerm
phi3t = var "phi3t" "\\Phi_{3t}" $ (n3 + (1-n3)**2*log(1-n3))/(36*pi * n3**2 * (1-n3)**2)*tensorThirdTerm

kTphi1, kTphi2, kTphi3, kTphi3t :: Expression Scalar
kTphi1 = var "kTphi1" "kT\\Phi_1" (integrate (kT*phi1))
kTphi2 = var "kTphi2" "kT\\Phi_2" (integrate (kT*phi2))
kTphi3 = var "kTphi3" "kT\\Phi_3" (integrate (kT*phi3))
kTphi3t = var "kTphi3t" "kT\\Phi_3" (integrate (kT*phi3t))

whitebear :: Expression Scalar
whitebear = var "whitebear" "F_{HS}" (kTphi1 + kTphi2 + kTphi3)

tensorwhitebear :: Expression Scalar
tensorwhitebear = var "tensorwhitebear" "F_{HSt}" (kTphi1 + kTphi2 + kTphi3t) -- $ integrate (kT*(phi1+phi2+phi3t))

phitot :: Expression RealSpace
phitot = var "phitot" "\\Phi" $ phi1 + phi2 + phi3

gSigmaS :: Expression RealSpace
gSigmaS = var "gSigmaS" "g_{\\sigma}^{S}" $
                   gSigmaS_helper phitot

gSigmaS_m2 :: Expression RealSpace
gSigmaS_m2 = var "gSigmaS" "g_{\\sigma,m2}^{S}" $
                      gSigmaS_helper phitot_m2

gSigmaS_helper :: Expression RealSpace -> Expression RealSpace
gSigmaS_helper phit = dAdR/(kT * n0**2 * 4*pi* (2*rad)**2)
    where dAdR = var "dAdR" "\\frac{dA}{dR}" $
                 kT*((d n3)*n2
                     - (d n2)*n2p'
                     - (d n1)*( n2p'/(4*pi*rad) + n0)
                     - (d n0) * ( n2p'/(4*pi*rad**2) + 2*n0/rad )
                     - (dv n2v `dot` n2vp' + dv n1v `dot` ( n2vp' - n2v/.rad)/(4*pi*rad)))
          d a = derive a 1 phit
          dv a = realspaceGradient a phit
          n2p' = -2*n2/rad
          n2vp' = (-2/rad).*n2v
          
nA :: Expression RealSpace
nA = "nA" === double_shell n / (4*pi*(2*rad)**2)

double_shell :: Expression RealSpace -> Expression RealSpace
double_shell x = ifft ( deltak * fft x)
  where deltak = var "delta2k" "\\delta_2(k)" $ smear * (4*pi) * (2*rad) * sin (2*kR) / k

gSigmaA :: Expression RealSpace
gSigmaA = gSigmaA_automagic

gSigmaA_automagic :: Expression RealSpace
gSigmaA_automagic = var "gSigmaA" "g_{\\sigma}^{A}" $
                    gSigmaA_helper phitot

gSigmaA_m2 :: Expression RealSpace
gSigmaA_m2 = var "gSigmaA" "g_{\\sigma,m2}^{A}" $
                      gSigmaA_helper phitot_m2

gSigmaA_by_hand :: Expression RealSpace
gSigmaA_by_hand = dAdR/(kT * n*shell_diam n )
    where dAdR = var "dAdR" "\\frac{dA}{dR}" $ factorize $
                 kT*n*( shell (dphi_dn3 - dphi_dn1/(4*pi*rad**2) - dphi_dn0/(2*pi*rad**3))
                        - shellPrime (dphi_dn2 + dphi_dn1/(4*pi*rad) + dphi_dn0/(4*pi*rad**2))
                        + vshellPrimedot (dphi_dn2v + 1/(4*pi*rad) .* dphi_dn1v)
                        + vshelldot (1/(4*pi*rad**2) .* dphi_dn1v) )
          dphi_dn0 = -log(1-n3)
          dphi_dn1 = n2/(1-n3)
          dphi_dn2 = n1/(1-n3) + (3*n2**2 - 3*sqr_n2v)*(n3+(1-n3)**2*log(1-n3))/(36*pi*n3**2*(1-n3)**2)
          dphi_dn3 = n0/(1-n3) + (n1*n2- n1v_dot_n2v)/(1-n3)**2 +
                     vectorThirdTerm*(1 - 2*(1-n3)*log(1-n3) - (1-n3))/(36*pi*n3**2*(1-n3)**2) -
                     vectorThirdTerm*(n3+(1-n3)**2*log(1-n3))/(36*pi*n3**2*(1-n3)**2)**2*(36*pi)*(2*n3*(1-n3)**2 - 2*n3**2*(1-n3))
          dphi_dn1v = -1/(1-n3) .* n2v
          dphi_dn2v = -1/(1-n3) .* n1v - 6*n2*(n3+(1-n3)**2*log(1-n3))/(36*pi*n3**2*(1-n3)**2) .* n2v

gSigmaA_helper :: Expression RealSpace -> Expression RealSpace
gSigmaA_helper phit = dAdR/(kT * n*shell_diam n )
    where dAdR = -- var "dAdR" "\\frac{dA}{dR}" $ 
                 factorize $
                 kT*n*( shell (d n3)
                        - shellPrime (d n2)
                        - ( shellPrime (d n1) + (shell (d n1))/rad ) / (4*pi*rad)
                        - ( shellPrime (d n0) + 2*(shell (d n0))/rad ) / (4*pi*rad**2)
                        + ( vshellPrimedot (dv n2v)  +
                            ( vshellPrimedot (dv n1v) + (vshelldot (dv n1v))/rad )/(4*pi*rad) ))
          d a = derive a 1 phit
          dv a = realspaceGradient a phit

correlation_gross :: Expression RealSpace
correlation_gross = (1 - 0.5*eta)/(1 - eta)**3
    where eta = step_diam n/8

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
