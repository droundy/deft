module SW_liquid
       (sw_liquid_n)
       where

import Expression
import IdealGas ( idealgas )
import WhiteBear ( gSigmaA, whitebear )

sigma :: Type a => Expression a
sigma = s_tex "sigma" "\\sigma"

epsilon :: Type a => Expression a
epsilon = s_tex "epsilon" "\\epsilon"

ksig :: Expression KSpace
ksig = k*sigma

lamksig :: Expression KSpace
lamksig = lambda*k*sigma

lambda :: Type a => Expression a
lambda = s_tex "lambda" "\\lambda"

eps4piok :: Expression KSpace
eps4piok = epsilon*4*pi/k

n :: Expression RealSpace
n = r_var "n"

gsigma :: Expression RealSpace
gsigma = substitute ("n" === r_var "x") n gSigmaA

sw_liquid_n :: Expression Scalar
sw_liquid_n = "ESW" === (substitute ("n" === r_var "x") (r_var "n") $
                         sw + whitebear + idealgas + 
                         ("external" === integrate (n * (r_var "Vext" - s_var "mu"))))

sw :: Expression Scalar
sw = var "sw" "F_{\\text{sw}}" $ epsilon * (0*sw0 + 0*sw1 + 0*sw2 + 0*sw3 + 0*sw4)
  where sw0 = var "sw0" "\\Sigma_0" $ integrate $ 0.5*n*n_g_phi0
        sw1 = var "sw1" "\\Sigma_1" $ integrate $ 0.5*n*n_g_phi1
        sw2 = var "sw2" "\\Sigma_2" $ integrate $ 0.5*n*n_g_phi2
        sw3 = var "sw3" "\\Sigma_3" $ integrate $ 0.5*n*n_g_phi3
        sw4 = var "sw4" "\\Sigma_4" $ integrate $ 0.5*n*n_g_phi4

n_g_phi0, n_g_phi1, n_g_phi2, n_g_phi3, n_g_phi4 :: Expression RealSpace
n_g_phi0 = "ngphi0" === convolve_xi0phi_with (n * gsigma)
n_g_phi1 = "ngphi1" === convolve_xi1phi_with (n * (k11*(gsigma-1) + k21*(gsigma-1)**2 +
                                  k31*(gsigma-1)**3 + k41*(gsigma-1)**4))
  where k11 = -1.754
        k21 = 0.027
        k31 = 0.838
        k41 = -0.178
n_g_phi2 = "ngphi2" === convolve_xi2phi_with (n * (k12*(gsigma-1) + k22*(gsigma-1)**2 +
                                  k32*(gsigma-1)**3 + k42*(gsigma-1)**4))
  where k12 = -2.243
        k22 = 4.403
        k32 = -2.48
        k42 = 0.363
n_g_phi3 = "ngphi3" === convolve_xi3phi_with (n * (k13*(gsigma-1) + k23*(gsigma-1)**2 +
                                  k33*(gsigma-1)**3 + k43*(gsigma-1)**4))
  where k13 = 0.207
        k23 = 0.712
        k33 = -1.952
        k43 = 1.046
n_g_phi4 = "ngphi4" === convolve_xi4phi_with (n * (k14*(gsigma-1) + k24*(gsigma-1)**2 +
                                  k34*(gsigma-1)**3 + k44*(gsigma-1)**4))
  where k14 = -0.002
        k24 = -0.164
        k34 = 0.324
        k44 = -0.162
        
convolve_xi0phi_with :: Expression RealSpace -> Expression RealSpace
convolve_xi0phi_with x = ifft ( xi0phik * fft x)
  where xi0phik = var "xi1phik" "\\tilde{\\xi_0\\phi}(k)" $
               -eps4piok * (sin(lamksig)-lamksig*cos(lamksig)-sin(ksig)+ksig*cos(ksig))/k**2

convolve_xi1phi_with :: Expression RealSpace -> Expression RealSpace
convolve_xi1phi_with x = ifft ( xi1phik * fft x)
  where xi1phik = var "xi1phik" "\\tilde{\\xi_1\\phi}(k)" $
               -eps4piok * ((2-lamksig*ksig*(lambda-1))*cos(lamksig)
                            -ksig*(sin(ksig)+(1-2*lambda)*sin(lamksig))-2*cos(ksig))/(k**2*ksig)

convolve_xi2phi_with :: Expression RealSpace -> Expression RealSpace
convolve_xi2phi_with x = ifft ( xi2phik * fft x)
  where xi2phik = var "xi2phik" "\\tilde{\\xi_2\\phi}(k)" $
               -eps4piok * (6*sin(ksig) 
                            + (ksig**2*(1-4*lambda+3*lambda**2)-6)*sin(lamksig)
                            - ksig*(4+lambda*(ksig**2*(lambda-1)**2-6))*cos(lamksig)
                            - 2*ksig*cos(ksig))/(k**2*ksig**2)

convolve_xi3phi_with :: Expression RealSpace -> Expression RealSpace
convolve_xi3phi_with x = ifft ( xi3phik * fft x)
  where xi3phik = var "xi3phik" "\\tilde{\\xi_3\\phi}(k)" $
               -eps4piok*(ksig*(6*sin(ksig) 
                                + (18-24*lambda+ksig**2*(lambda-1)**2*(4*lambda-1))*sin(lamksig))
                          + (6*ksig**2*(lambda-1)*(2*lambda-1)
                             - lamksig*ksig**3*(lambda-1)**3-24)*cos(lamksig)
                          + 24*cos(ksig))/(k**2*ksig**3)

convolve_xi4phi_with :: Expression RealSpace -> Expression RealSpace
convolve_xi4phi_with x = ifft ( xi4phik * fft x)
  where xi4phik = var "xi4phik" "\\tilde{\\xi_4\\phi}(k)" $
               -eps4piok *(lambda*(36*sin lamksig/(k**2*ksig**2)
                                   - 60*lambda*sin lamksig/(k**2*ksig**2)
                                   + (lambda-1)**2*(5*lambda-1)*sin lamksig/k**2)
                           - (36*sin lamksig/(k**2*ksig**2)
                              - 60*lambda*sin lamksig/(k**2*ksig**2)
                              + (lambda-1)**2*(5*lambda-1)*sin lamksig/k**2)
                           + (120*sin lamksig - 120*sin ksig)/(k**2*ksig**4)
                           + 24*cos ksig/(k**2*ksig**3)
                           - 24*(5*lambda-4)*cos lamksig/(k**2*ksig**3)
                           + 4*(lambda-1)**2*(5*lambda-2)*cos lamksig/(k**2*ksig)
                           - lamksig*(lambda-1)**4*cos lamksig/(k**2))
