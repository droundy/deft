-- This is the SAFT functional defined in the Krebs and Roundy
-- paper (to be published), and described in papers/water-SAFT.  It is
-- in its own separate module, so that we will leave it unmodified for
-- posterity (and for future comparisons).

-- This module also defines some functions from Yu and Wu (2002),
-- which are used in the Hughes functional.

module WaterSaft
       ( eta_for_dispersion, lambda_dispersion, a1, a2, eta_effective,
         saft_dispersion, saft_association, water_saft, water_entropy,
         water_X, mu, water_saft_by_hand,
         water_saft_n, water_saft_by_hand_n, homogeneous_water_saft_n, homogeneous_water_saft_by_hand_n)
       where

import FMT ( rad, n )
import WhiteBear ( whitebear, kT, gSigmaA, gSigmaA_by_hand, nA )
import IdealGas ( idealgas )
import Expression

water_saft_n, water_saft_by_hand_n, homogeneous_water_saft_n, homogeneous_water_saft_by_hand_n :: Expression Scalar
water_saft_n = substitute n (r_var "n") water_saft
water_saft_by_hand_n = substitute n (r_var "n") water_saft_by_hand

homogeneous_water_saft_n = makeHomogeneous water_saft_n
homogeneous_water_saft_by_hand_n = makeHomogeneous water_saft_by_hand_n

mu :: Type a => Expression a
mu = s_tex "mu" "\\mu"

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

saft_association = "Fassoc" === integrate (4*kT*n*(log water_X - water_X/2 + 1/2))

kappa_association, epsilon_association :: Type a => Expression a
kappa_association = s_tex "kappa_association" "\\kappa_a"
epsilon_association = s_tex "epsilon_association" "\\epsilon_a"

water_X, deltasaft, a1, a2 :: Expression RealSpace
water_X = "X" === (sqrt(1 + 8*nA*deltasaft) - 1) / (4*nA*deltasaft)

deltasaft = var "deltasaft" "{\\Delta}" $ gSW*kappa_association*fmayer
  where fmayer = "fmayer" === exp(epsilon_association/kT)-1

gSW :: Expression RealSpace
gSW = "gSW" ===
      gSigmaA
      + (1/4/kT)*(da1_detad - lambda_dispersion/(3*eta_for_dispersion)*da1_dlambda)

water_saft :: Expression Scalar
water_saft = "FSAFT" === (idealgas + whitebear + saft_association + saft_dispersion + integrate (n*mu))

water_entropy :: Expression Scalar
water_entropy = "SSAFT" === (d "Sig" idealgas + d "Shs" whitebear + d "Sass" saft_association + d "Sdisp" saft_dispersion)
  where d nn f = nn === - scalarderive (ES kT) f




saft_association_by_hand :: Expression Scalar
saft_association_by_hand = "Fassoc" === integrate (4*kT*n*(log water_X_by_hand - water_X_by_hand/2 + 1/2))

water_X_by_hand, deltasaft_by_hand :: Expression RealSpace
water_X_by_hand = "X" === (sqrt(1 + 8*nA*deltasaft_by_hand) - 1) / (4*nA*deltasaft_by_hand)
deltasaft_by_hand = var "deltasaft" "{\\Delta}" $ gSW_by_hand*kappa_association*boltz
  where boltz = "boltz" === exp(epsilon_association/kT)-1

gSW_by_hand :: Expression RealSpace
gSW_by_hand = "gSW" ===
      gSigmaA_by_hand
      + (1/4/kT)*(da1_detad - lambda_dispersion/(3*eta_for_dispersion)*da1_dlambda)


water_saft_by_hand :: Expression Scalar
water_saft_by_hand = "FSAFT" === (idealgas + whitebear + saft_association_by_hand + saft_dispersion + integrate (n*mu))
