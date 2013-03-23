-- This is the SAFT functional defined in the Hughes, Krebs and Roundy
-- paper (to be published), and described in papers/water-SAFT.  It is
-- in its own separate module, so that we will leave it unmodified for
-- posterity (and for future comparisons).

-- This module also defines some functions from Yu and Wu (2002),
-- which are used in the Hughes functional.

module HughesSaft
       ( yuwu_zeta, yuwu_correlation,
         eta_for_dispersion, lambda_dispersion, a1, a2, eta_effective,
         saft_dispersion, saft_association, saft_fluid, saft_entropy, mu )
       where

import FMT ( rad, n, n0, n2, n3, sqr_n2v )
import WhiteBear ( whitebear, kT )
import IdealGas ( idealgas )
import Expression

mu :: Type a => Expression a
mu = s_tex "mu" "\\mu"

yuwu_zeta :: Expression RealSpace
yuwu_zeta = var "zeta_yuwu" "{\\zeta}" $ (1 - sqr_n2v/n2**2)

yuwu_correlation :: Expression RealSpace
yuwu_correlation = var "ghsyuwu" "g_{HS}^{\\textit{YuWu}}" ghs
    where ghs = 1/(1-n3) + rad/2*n2*yuwu_zeta/(1-n3)**2 +
                rad**2/18*n2**2*yuwu_zeta/(1-n3)**3
    --where ghs = (1 + 0.5*(rad*n2/3/(1-n3))*yuwu_zeta*(3 + rad*n2/3/(1-n3)))/(1-n3)

lambda_dispersion, epsilon_dispersion :: Type a => Expression a
lambda_dispersion = s_tex "lambda_dispersion" "\\lambda_d"
epsilon_dispersion = s_tex "epsilon_dispersion" "\\epsilon_d"

eta_for_dispersion, eta_effective :: Expression RealSpace
--eta_for_dispersion = var "eta_d" "{\\eta_d}" $ (4*pi*rad**3/3)*n
eta_for_dispersion = var "eta_d" "{\\eta_d}" $
                     ((4*pi*rad**3/3)*ifft (exp (-0.5*k**2*(2*length_scaling*lambda_dispersion*rad)**2) * fft n))
                       where length_scaling = s_tex "length_scaling" "l_s"

-- The following equation is equation 36 in Gil-Villegas 1997 paper.
eta_effective = var "eta_eff" "\\eta_{\\textit{eff}}" $
                (c1 + c2*eta_for_dispersion + c3*eta_for_dispersion**2)*eta_for_dispersion
  where c1 = "c1" === 2.25855 - 1.50349*lambda_dispersion + 0.249434*lambda_dispersion**2
        c2 = "c2" === -0.669270 + 1.40049*lambda_dispersion - 0.827739*lambda_dispersion**2
        c3 = "c3" === 10.1576 - 15.0427*lambda_dispersion + 5.30827*lambda_dispersion**2

saft_dispersion, saft_association :: Expression Scalar
saft_dispersion = "Fdisp" === (("a1integrated" === integrate (n*a1)) + ("a2integrated" === integrate (n*a2/kT)))

a1 = "a1" === -4*(lambda_dispersion**3-1)*epsilon_dispersion*eta_for_dispersion*ghs
  where ghs = var "ghs" "g_{HS}" $ (1 - eta_effective/2)/(1-eta_effective)**3

da1_dlambda, da1_detad :: Expression RealSpace
da1_dlambda = substitute (r_var "etad") eta_for_dispersion $ derive (lambda_dispersion :: Expression RealSpace) (1 :: Expression RealSpace) $ substitute eta_for_dispersion (r_var "etad") a1

da1_detad = derive eta_for_dispersion 1 a1

a2 = "a2" === 0.5*khs*epsilon_dispersion*eta_for_dispersion*da1_detad
     where khs = var "KHS" "{\\kappa_{HS}}" $ (1 - eta_for_dispersion)**4/(1 + 4*(eta_for_dispersion + eta_for_dispersion**2))

saft_association = "Fassoc" === integrate (4*kT*n0*yuwu_zeta*(log xsaft - xsaft/2 + 1/2))

kappa_association, epsilon_association :: Type a => Expression a
kappa_association = s_tex "kappa_association" "\\kappa_a"
epsilon_association = s_tex "epsilon_association" "\\epsilon_a"

xsaft, deltasaft, a1, a2 :: Expression RealSpace
xsaft = "X" === (sqrt(1 + 8*n0*yuwu_zeta*deltasaft) - 1) / (4*n0*yuwu_zeta*deltasaft)

deltasaft = var "deltasaft" "{\\Delta}" $ gSW*kappa_association*boltz
  where boltz = "boltz" === exp(epsilon_association/kT)-1

gSW :: Expression RealSpace
gSW = "gSW" ===
      yuwu_correlation
      + (1/4/kT)*(da1_detad - lambda_dispersion/(3*eta_for_dispersion)*da1_dlambda)

saft_fluid :: Expression Scalar
saft_fluid = "FSAFT" === (idealgas + whitebear + saft_association + saft_dispersion + integrate (n*mu))

saft_entropy :: Expression Scalar
saft_entropy = "SSAFT" === (d "Sig" idealgas + d "Shs" whitebear + d "Sass" saft_association + d "Sdisp" saft_dispersion)
  where d nn f = nn === - scalarderive (ES kT) f
